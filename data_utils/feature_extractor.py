MAX_LEN=30
import os
import pickle
import re
import numpy as np
from collections import defaultdict
import argparse

class data_preprocess(object):
  """
  This is a helper class designed for multi-process get_spectrum
  """

  def dict_to_pkl(self, fname, dict):

      with open(fname, mode='wb') as handle:
          pickle.dump(dict, handle, pickle.HIGHEST_PROTOCOL)

  def __init__(self, input_spectrum_file, input_feature_file):
    self.MZ_MAX = 3000
    self.SPECTRUM_RESOLUTION = 50
    self.MZ_SIZE = MZ_SIZE = int(self.MZ_MAX * self.SPECTRUM_RESOLUTION)
    self.neighbor_size = 5

    self.dia_window = 10.0 # the window size of MS2 scan in Dalton

    self.input_spectrum_file = input_spectrum_file
    self.input_feature_file = input_feature_file
    self.output_file = ''

    # split data into batches
    self.feature_index_list =[]


    ### store file location of each feature for random access
    self.feature_location_list = []

    # store the file location of all spectra for random access
    self.spectrum_location_dict =  {}
    self.spectrum_rtinseconds_dict = {}

    # record the status of spectra that have been read
    self.feature_count = {"total": 0,
                          "read": 0,
                          "skipped": 0,
                          "skipped_mass": 0}
    self.spectrum_count = 0
    self.precursor_info = defaultdict(list)
    self.col_feature_id = 0
    self.col_precursor_mz = 1
    self.col_precursor_charge = 2
    self.col_rt_mean = 3
    self.col_raw_sequence = 4
    self.col_scan_list = 5
    self.col_ms1_list = 6
    self.col_feature_area = 7
    self.col_num = 8
    # predicted file column format
    self.pcol_feature_id = 0
    self.pcol_feature_area = 1
    self.pcol_sequence = 2
    self.pcol_score = 3
    self.pcol_position_score = 4
    self.pcol_precursor_mz = 5
    self.pcol_precursor_charge = 6
    self.pcol_protein_id = 7
    self.pcol_scan_list_middle = 8
    self.pcol_scan_list_original = 9
    self.pcol_score_max = 10
    self.mass_H = 1.0078
    self.mass_H2O = 18.0106
    self.mass_NH3 = 17.0265
    self.mass_N_terminus = 1.0078
    self.mass_C_terminus = 17.0027
    self.mass_CO = 27.9949

  def get_location(self):

    print("".join(["="] * 80)) # section-separating line
    print("Feature Extraction: get_location()")

    ### store file location of each spectrum for random access {scan:location}
    ### since mgf file can be rather big, cache the locations for each spectrum mgf file.
    spectrum_location_file = self.input_spectrum_file + '.locations.pkl'
    if os.path.exists(spectrum_location_file):
      print("Feature Extraction: read cached spectrum locations")
      with open(spectrum_location_file, 'rb') as fr:
        data = pickle.load(fr)
        self.spectrum_location_dict, self.spectrum_rtinseconds_dict, self.spectrum_count = data
    else:
      print("Feature Extraction: build spectrum location from scratch")
      spectrum_location_dict = {}
      spectrum_rtinseconds_dict = {}
      line = True
      while line:
        current_location = self.input_spectrum_handle.tell()
        line = self.input_spectrum_handle.readline()
        if "BEGIN IONS" in line:
          spectrum_location = current_location
        elif "SCANS=" in line:
          scan = re.split('=|\r\n', line)[1]
          spectrum_location_dict[scan] = spectrum_location
        elif "RTINSECONDS=" in line:
          rtinseconds = float(re.split('=|\r\n', line)[1])
          spectrum_rtinseconds_dict[scan] = rtinseconds
      self.spectrum_location_dict = spectrum_location_dict
      self.spectrum_rtinseconds_dict = spectrum_rtinseconds_dict
      self.spectrum_count = len(spectrum_location_dict)
      with open(spectrum_location_file, 'wb') as fw:
        pickle.dump((self.spectrum_location_dict, self.spectrum_rtinseconds_dict, self.spectrum_count), fw)

    ### store file location of each feature for random access
    # skip header line
    _ = self.input_feature_handle.readline()
    line = True
    while line:
      feature_location = self.input_feature_handle.tell()
      self.feature_location_list.append(feature_location)
      line = self.input_feature_handle.readline()
    self.feature_location_list = self.feature_location_list[:-1]
    self.feature_location_list = self.feature_location_list
    self.feature_count["total"] = len(self.feature_location_list)
    self.feature_index_list = range(self.feature_count["total"])

    print("spectrum_count = {0:d}".format(self.spectrum_count))
    print("feature_count[total] = {0:d}".format(self.feature_count["total"]))

  def open_input(self):

    print("".join(["="] * 80)) # section-separating line
    print("Feature Extraction: open_input()")

    self.input_spectrum_handle = open(self.input_spectrum_file, 'r')
    self.input_feature_handle = open(self.input_feature_file, 'r')

  def _print_prediction_header(self):
    print("".join(["="] * 80)) # section-separating line
    print("WorkerIO: _print_prediction_header()")

    header_list = ["feature_id",
                   "feature_area",
                   "predicted_sequence",
                   "predicted_score",
                   "predicted_position_score",
                   "precursor_mz",
                   "precursor_charge",
                   "protein_access_id",
                   "scan_list_middle",
                   "scan_list_original",
                   "predicted_score_max"]
    header_row = "\t".join(header_list)
    print(header_row, file=self.output_handle, end="\n")
  def open_output(self):

    print("".join(["="] * 80)) # section-separating line
    print("Feature Extraction: open_output()")

    self.output_handle = open(self.output_file, 'w')
    self._print_prediction_header()


  def get_spectrum(self, feature_index_batch, input_feature_file_handle, input_spectrum_file_handle):
    spectrum_list = []
    for feature_index in feature_index_batch:
      # parse a feature
      feature_location = self.feature_location_list[feature_index]
      feature_id, feature_area, precursor_mz, precursor_charge, rt_mean, raw_sequence, scan_list, ms1_list = self._parse_feature(feature_location, input_feature_file_handle)
      # skip if precursor_mass > MZ_MAX
      precursor_mass = precursor_mz * precursor_charge - self.mass_H * precursor_charge
      if precursor_mass > self.MZ_MAX:
        continue

      # parse and process spectrum
      (neighbor_right_count,
       neighbor_size_half,
       scan_list_middle,
        mz_list,
       intensity_list,
       ms1_profile) = self._parse_spectrum(precursor_mz, precursor_mass, rt_mean, scan_list, ms1_list, input_spectrum_file_handle)
      # update dataset
      spectrum = {
                  "raw_sequence": raw_sequence,
                  "precursor_mass": precursor_mass,
                  "precursor_mz": precursor_mz,
                  "precursor_charge": precursor_charge,
                  "ms1_profile": ms1_profile,
                  "mz_list":mz_list,
                  "intensity_list": intensity_list,
                  "scan_list_middle": scan_list_middle,
                  "neighbor_right_count":neighbor_right_count,
                  "neighbor_size_half":neighbor_size_half}
      spectrum_list.append(spectrum)

    return spectrum_list

  def _parse_feature(self, feature_location, input_file_handle):
    input_file_handle.seek(feature_location)
    line = input_file_handle.readline()
    line = re.split(',|\r|\n', line)
    feature_id = line[self.col_feature_id]
    feature_area = 0#float(line[deepnovo_config.col_feature_area])
    precursor_mz = float(line[self.col_precursor_mz])
    precursor_charge = float(line[self.col_precursor_charge])
    rt_mean = float(line[self.col_rt_mean])
    raw_sequence = line[self.col_raw_sequence]
    scan_list = re.split(';', line[self.col_scan_list])
    ms1_list = re.split(';', line[self.col_ms1_list])
    assert len(scan_list) == len(ms1_list), "Error: scan_list and ms1_list not matched."

    return feature_id, feature_area, precursor_mz, precursor_charge, rt_mean, raw_sequence, scan_list, ms1_list

  def _parse_spectrum_ion(self, input_file_handle):
    mz_list = []
    intensity_list = []
    line = input_file_handle.readline()
    while not "END IONS" in line:
      mz, intensity = re.split(' |\n', line)[:2]
      mz_float = float(mz)
      intensity_float = float(intensity)
      # skip an ion if its mass > MZ_MAX
      if mz_float > self.MZ_MAX:
        line = input_file_handle.readline()
        continue
      mz_list.append(mz_float)
      intensity_list.append(intensity_float)
      line = input_file_handle.readline()

    return mz_list, intensity_list

  def _parse_spectrum(self, precursor_mz, precursor_mass, rt_mean, scan_list, ms1_list, input_file_handle):
    mz_five_list = []
    int_five_list = []
    neighbor_count = len(scan_list)
    best_scan_index = None
    best_distance = float('inf')
    for scan_index, scan in enumerate(scan_list):
      try:
        distance = abs(self.spectrum_rtinseconds_dict[scan] - rt_mean)
      except:
        distance = abs(self.spectrum_rtinseconds_dict[scan+ '\n'] - rt_mean)
      if distance < best_distance:
        best_distance = distance
        best_scan_index = scan_index
    neighbor_center = best_scan_index
    neighbor_left_count = neighbor_center
    neighbor_right_count = neighbor_count - neighbor_left_count - 1
    neighbor_size_half = self.neighbor_size // 2
    neighbor_left_count = min(neighbor_left_count, neighbor_size_half)
    neighbor_right_count = min(neighbor_right_count, neighbor_size_half)
    ### parse and add neighbor spectra
    scan_list_middle = []
    ms1_intensity_list_middle = []
    for index in range(neighbor_center - neighbor_left_count, neighbor_center + neighbor_right_count + 1):
      scan = scan_list[index]
      scan_list_middle.append(scan)
      ms1_entry = ms1_list[index]
      ms1_intensity = float(re.split(':', ms1_entry)[1])
      ms1_intensity_list_middle.append(ms1_intensity)
    ms1_intensity_max = max(ms1_intensity_list_middle)
    assert ms1_intensity_max > 0.0, "Error: Zero ms1_intensity_max"
    ms1_intensity_list_middle = [x/ms1_intensity_max for x in ms1_intensity_list_middle]
    for scan, ms1_intensity in zip(scan_list_middle, ms1_intensity_list_middle):
      try:
          spectrum_location = self.spectrum_location_dict[scan]
      except:
          spectrum_location = self.spectrum_location_dict[scan+'\n']
      input_file_handle.seek(spectrum_location)
      # parse header lines
      line = input_file_handle.readline()
      assert "BEGIN IONS" in line, "Error: wrong input BEGIN IONS"
      line = input_file_handle.readline()
      assert "TITLE=" in line, "Error: wrong input TITLE="
      line = input_file_handle.readline()
      assert "PEPMASS=" in line, "Error: wrong input PEPMASS="
      line = input_file_handle.readline()
      assert "CHARGE=" in line, "Error: wrong input CHARGE="
      line = input_file_handle.readline()
      assert "SCANS=" in line, "Error: wrong input SCANS="
      line = input_file_handle.readline()
      assert "RTINSECONDS=" in line, "Error: wrong input RTINSECONDS="
      # parse fragment ions
      mz_list, intensity_list = self._parse_spectrum_ion(input_file_handle)
      mz_five_list.append(mz_list)
      int_five_list.append(intensity_list)

    # ms1_profile
    for x in range(neighbor_size_half - neighbor_left_count):
      ms1_intensity_list_middle = [0.0] + ms1_intensity_list_middle
    for x in range(neighbor_size_half - neighbor_right_count):
      ms1_intensity_list_middle = ms1_intensity_list_middle + [0.0]
    assert len(ms1_intensity_list_middle) == self.neighbor_size, "Error: ms1 profile"
    ms1_profile = np.array(ms1_intensity_list_middle)

    return neighbor_right_count, neighbor_size_half, scan_list_middle, mz_five_list, int_five_list, ms1_profile

  def read_precursor_info(self, feature_index, feature_fr, spectrum_fr):

      spectrum_list = self.get_spectrum([feature_index], feature_fr, spectrum_fr)
      if not spectrum_list:
          return None

      spectrum = spectrum_list[0]
      raw_sequence = spectrum["raw_sequence"]
      spec_info = {
          "precursor_mz": spectrum["precursor_mz"],
          "precursor_charge": spectrum["precursor_charge"],
          "scan_list_middle": spectrum["scan_list_middle"],
          "ms1": spectrum["ms1_profile"],
          "mz_list": spectrum["mz_list"],
          "int_list": spectrum["intensity_list"],
          "neighbor_right_count": spectrum["neighbor_right_count"],
          "neighbor_size_half": spectrum["neighbor_size_half"]}
      ptms = ['C(+57.02)']
      ptm_casanovo_format = 'C+57.021'
      for ptm in ptms:

          if ptm in raw_sequence:
              raw_sequence = raw_sequence.replace(ptm, ptm_casanovo_format)
      self.precursor_info[raw_sequence].append(spec_info)
      return
  def get_precursor_info(self):
      with open(self.input_feature_file, 'r') as feature_fr:
          with open(self.input_spectrum_file, 'r') as spectrum_fr:
              for feature_index in self.feature_index_list:
                  self.read_precursor_info(feature_index, feature_fr, spectrum_fr)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process spectrum and feature files to generate a feature file in pickle format for Transformer_DIA or DiffNovo_DIA.")
    parser.add_argument("--feature_file", required=True, help="Path to the feature file (.csv)")
    parser.add_argument("--spectrum_file", required=True, help="Path to the spectrum file (.mgf)")
    parser.add_argument("--output_file", required=True, help="Path to save the output pickle file (.pkl)")

    args = parser.parse_args()

    # Validate input files
    if not os.path.isfile(args.feature_file) or not args.feature_file.endswith(".csv"):
        raise ValueError("Feature file must exist and have a .csv extension.")
    if not os.path.isfile(args.spectrum_file) or not args.spectrum_file.endswith(".mgf"):
        raise ValueError("Spectrum file must exist and have a .mgf extension.")

    print(
        "This script processes spectrum and feature files to generate a feature file in pickle format, "
        "required for the Transformer_DIA or DiffNovo_DIA model.\n"
        "Generated features:\n"
        "- Keys: Peptide sequences\n"
        "- Values: [precursor_mz, precursor_charge, scan_list_middle, ms1, mz_list, int_list, neighbor_right_count, neighbor_size_half]"
    )

    # Run preprocessing
    feature_extraction = data_preprocess(args.spectrum_file, args.feature_file)
    feature_extraction.open_input()
    feature_extraction.get_location()
    feature_extraction.get_precursor_info()
    feature_extraction.dict_to_pkl(args.output_file, feature_extraction.precursor_info)

