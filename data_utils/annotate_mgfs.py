import pandas as pd
import re
import os

from collections import defaultdict
import pickle
import numpy as np
import argparse
import sys

MZ_MAX = 3000.0
mass_H = mass_H = 1.0078

def normalize_int(int_list):
    return np.sqrt(np.sqrt(int_list))

def get_scan_int(scan_set, spectrum_location_dict, input_spectrum_handle):
    top_five_scans = []
    #for each scan number there is only one spectrum in mgf file
    # but there are multiple peptide sequences in the feature file for each specturm
    scan_int = {}
    for scan  in scan_set:
        spectrum_location = spectrum_location_dict[scan]
        input_spectrum_handle.seek(spectrum_location)
        # parse header lines
        line = input_spectrum_handle.readline()
        assert "BEGIN IONS" in line, "Error: wrong input BEGIN IONS"
        line = input_spectrum_handle.readline()
        assert "TITLE=" in line, "Error: wrong input TITLE="
        line = input_spectrum_handle.readline()
        assert "PEPMASS=" in line, "Error: wrong input PEPMASS="
        line = input_spectrum_handle.readline()
        assert "CHARGE=" in line, "Error: wrong input CHARGE="
        line = input_spectrum_handle.readline()
        assert "SCANS=" in line, "Error: wrong input SCANS="
        line = input_spectrum_handle.readline()
        assert "RTINSECONDS=" in line, "Error: wrong input RTINSECONDS="
        # parse fragment ions
        mz_list, intensity_list = _parse_spectrum_ion()
        normal_it = normalize_int(intensity_list)
        scan_int[scan] = np.sum(normal_it)
    return scan_int



def extract_scans_feature_f(feature_f):
    feature_df = pd.read_csv(feature_f, usecols=[0, 1, 2, 3, 4, 5])
    scans_set = set()
    num_scans = 0
    for _, row in feature_df.iterrows():
        spec_group_id, mz, z, rt_mean, seq, scans = row
        scan_list = scans.split(';')
        for scan in scan_list:
            scans_set.add(scan)

    return scans_set
def _parse_spectrum_ion(input_file_handle):
    """TODO(nh2tran): docstring."""

    #~ print("".join(["="] * 80)) # section-separating line
    #~ print("WorkerIO: _parse_spectrum_ion()")

    # ion
    mz_list = []
    intensity_list = []
    line = input_file_handle.readline()
    while not "END IONS" in line:
      mz, intensity = re.split(' |\n', line)[:2]
      mz_float = float(mz)
      intensity_float = float(intensity)
      # skip an ion if its mass > MZ_MAX
      if mz_float > MZ_MAX:
        line = input_file_handle.readline()
        continue
      mz_list.append(mz_float)
      intensity_list.append(intensity_float)
      line = input_file_handle.readline()

    return mz_list, intensity_list	

def annotate_mgf_one_spectrum_one_precursor (spect_f, spectrum_location_dict, precursor_dict, model, all_scans = False):
	        
	        out_f = spect_f.split('.mgf')[0]+'_one_scan_one_spec_all_scans.mgf' if all_scans else spect_f.split('.mgf')[0]+'_one_scan_one_spec_middle_scans_annotation.mgf'
	        with open (spect_f, mode='r') as input_file_handle:
                      with open(out_f, mode='w') as output_handle:
      			#'F3:82340': [{1162.658: (581.32884, 2.0, 111.0, 'TYLPAGQSVLL')}, {1305.716: (652.858, 2.0, 112.0, 'PTEVLAGTQLLY')}],
                            for scan, prec_list in precursor_dict.items():
                                masses = []
                                for prec in prec_list:
                                    #print(prec)
                                    masses.append(prec[0])
                                min_mass_dist = float('inf')
                                line_buffer = []
                                try:
                                   spectrum_location = spectrum_location_dict[scan]
                                except:
                                   spectrum_location = spectrum_location_dict[scan + '\n']
                                input_file_handle.seek(spectrum_location)
                                line = input_file_handle.readline()
                                line_buffer.append(line)
                                assert "BEGIN IONS" in line, "Error: wrong input BEGIN IONS"
                                line = input_file_handle.readline()
                                line_buffer.append(line)
                                assert "TITLE=" in line, "Error: wrong input TITLE="
                                line = input_file_handle.readline()
                                #print(line)
                                best_mass = 0
                                for mass in masses:
                                    if min_mass_dist > abs(float(line.replace("PEPMASS=", "").strip('\n')) - mass):
                                       best_mass = mass
                                print(f'best_mass={best_mass}, masses={masses}, spec_mass={line.replace("PEPMASS=", "")}')
                                pepmass= None
                                charge = None
                                rt = None
                                seq = None 
                                #prec = (mz, z, rt_mean, seq)
                                for prec in precursor_dict[scan]:
                                    #print(scan, prec)
                                    if best_mass == prec[0]:
                                            pepmass, charge, rt, seq = best_mass, prec[1], prec[2], prec[3]
                                #print(f'pepmass={pepmass}, charge={charge}, rt={rt}, seq={seq}')
                                # read the mass from mgf file
                                #line_buffer.append(line)
                                #read pepmass from feature file
                                #line_buffer.append('PEPMASS=' + str(round(precursor[1], 5)) + '\n')
                                line_buffer.append('PEPMASS=' + str(round(pepmass, 5)) + '\n')
                                #line_buffer.append(line)
                                assert "PEPMASS=" in line, "Error: wrong input PEPMASS="
                                line = input_file_handle.readline()
                                #if model == 'pepnet':
                                #   line_buffer.append('CHARGE=' + str(int(precursor[2])) + '\n')
                                #else:
                                #   line_buffer.append('CHARGE=' + str(precursor[2]) + '\n')
                                #line_buffer.append(line)
                                #line_buffer.append('CHARGE=' + str(int(precursor[2])) + '\n')
                                line_buffer.append('CHARGE=' + str(int(charge)) + '\n')
                                assert "CHARGE=" in line, "Error: wrong input CHARGE="
                                line = input_file_handle.readline()
                                #if line in precursor[0]: #print (scan, precursor[0]) # precursor[0] contains middle scan list
                                #if model == 'transdia':
                                #   print('model is transdia')
                                #middle_scan_label = ';'.join(middle_scan_list)
                                #line_buffer.append('SCANS=' + middle_scan_label + '\n')
                                #else:
                                line_buffer.append('SCANS=' + scan + '\n')
                                assert "SCANS=" in line, "Error: wrong input SCANS="
                                line = input_file_handle.readline()
                                #line_buffer.append(line)
                                #line_buffer.append('RTINSECONDS=' + str(round(precursor[3], 4)) + '\n')
                                line_buffer.append('RTINSECONDS=' + str(round(rt, 4)) + '\n')
                                assert "RTINSECONDS=" in line, "Error: wrong input RTINSECONDS="
                                #add the sequence int the precursor dictionary
                                #seq =  precursor[4]
                                line_buffer.append('SEQ=' + seq + '\n')
                                #if model is pepnet, we need to add seq to the title
                                if model == 'pepnet':
                                   line_buffer[1] = f"TITLE={seq}" + '\n'
                                # write the metadata in int ouput file
                                for l in line_buffer:
                                    output_handle.write(l)
                                while line and not ("END IONS" in line):
                                    line = input_file_handle.readline()
                                    output_handle.write(line)
                                output_handle.write("\n")





def annotate_mgf_one_spectrum_multiple_precursor(spect_f, spectrum_location_dict, precursor_dict, model='', out_name= '',all_scans = 'all'):

  #spectrum_location_dict:
  #is a dictionary of scans as keys and locations of each spectrum as value
  #{scan:locations},
  #precursors_dict:
  #is a dictionary of precursors information
  #{scan: list of precursors [(mz, z, rtime, seq)]}

  #out_f = spect_f.split('.mgf')[0]+'_one_pep_one_spec_all_scans.mgf' if all_scans else spect_f.split('.mgf')[0]+'one_pep_one_spec_middle_scans_annotation.mgf'
  print(f'annotating model ={model}, with policy of = {all_scans}')
  out_f = out_name  #spect_f.split('.mgf')[0]+'_'+model+ '_'+ all_scans +'.mgf'
  with open (spect_f, mode='r') as input_file_handle:
    with open(out_f, mode='w') as output_handle:
      for scan, precursor_list in precursor_dict.items():
        for precursor in precursor_list:
          #precursor [middel_scan_list, mz, z, rt_mean, seq]
          # [(['F1:2022', 'F1:2049', 'F1:2076', 'F1:2103', 'F1:2130'], 444.21304000000003, 2, 3.1655995999999997, 'SC+57.021HTAVGR')]
          # each precursor is a tuple of (middle_scan_list, mz, charge, rtime, seq)
          scan_list = precursor[0]
          #print(f'scan_list{scan_list}')
          mass_dist = 'info' 
          for scan in scan_list:
            line_buffer = []
            try:
              spectrum_location = spectrum_location_dict[scan]
            except:
              spectrum_location = spectrum_location_dict[scan + '\n']

            input_file_handle.seek(spectrum_location)
            # parse header lines
            line = input_file_handle.readline()
            line_buffer.append(line)
            assert "BEGIN IONS" in line, "Error: wrong input BEGIN IONS"

            line = input_file_handle.readline()
            line_buffer.append(line)
            assert "TITLE=" in line, "Error: wrong input TITLE="

            line = input_file_handle.readline()
            # read the mass from mgf file
            #line_buffer.append(line)
            #read pepmass from feature file
       
            line_buffer.append('PEPMASS=' + str(round(precursor[1], 5)) + '\n')
            #line_buffer.append(line)
            assert "PEPMASS=" in line, "Error: wrong input PEPMASS="

            line = input_file_handle.readline()
            #if model == 'pepnet':
             #   line_buffer.append('CHARGE=' + str(int(precursor[2])) + '\n')
            #else:
             #   line_buffer.append('CHARGE=' + str(precursor[2]) + '\n')

            #line_buffer.append(line)
            line_buffer.append('CHARGE=' + str(int(precursor[2])) + '\n')
            assert "CHARGE=" in line, "Error: wrong input CHARGE="

            line = input_file_handle.readline()
            #if line in precursor[0]: #print (scan, precursor[0]) # precursor[0] contains middle scan list
            #if model == 'transdia':
             #   print('model is transdia')
                #middle_scan_label = ';'.join(middle_scan_list)
                #line_buffer.append('SCANS=' + middle_scan_label + '\n')
            #else:
            line_buffer.append('SCANS=' + scan + '\n')
            assert "SCANS=" in line, "Error: wrong input SCANS="

            line = input_file_handle.readline()
            #line_buffer.append(line)
            line_buffer.append('RTINSECONDS=' + str(round(precursor[3], 4)) + '\n')
            assert "RTINSECONDS=" in line, "Error: wrong input RTINSECONDS="

            #add the sequence int the precursor dictionary
            seq =  precursor[4]
            if 'C(+57.02)' in seq:
                seq=seq.replace('C(+57.02)', 'C+57.021')
                #print(seq)
            line_buffer.append('SEQ=' + seq + '\n')
            #if model is pepnet, we need to add seq to the title
            if model == 'pepnet':
                line_buffer[1] = f"TITLE={seq}" + '\n'
            # write the metadata in int ouput file
            for l in line_buffer:
              output_handle.write(l)
            while line and not ("END IONS" in line):
              line = input_file_handle.readline()
              output_handle.write(line)
            output_handle.write("\n")

def add_seq_to_title_pepnet(path, spec_f):

  with open (path+spec_f, mode='r') as input_file_handle:
    #with open(path+'annotated_middle_scan_'+spec_f, mode='w') as output_handle:
    with open(path+'test_plasma_pepnet.mgf', mode='w') as output_handle:

      line = input_file_handle.readline()
      while line:

            # each precursor is a tuple of (spec_grouop_id, mz, charge, rtime, seq) like(583.28314, 2, 11.1451, 'AATGEC+57.021TATVGK')
            line_buffer = []
            assert "BEGIN IONS" in line, "Error: wrong input BEGIN IONS"
            line_buffer.append(line)
            line = input_file_handle.readline()
            assert "TITLE=" in line, "Error: wrong input TITLE="
            line_buffer.append(line)
            line = input_file_handle.readline()
            assert "PEPMASS=" in line, "Error: wrong input PEPMASS="
            line_buffer.append(line)
            line = input_file_handle.readline()
            assert "CHARGE=" in line, "Error: wrong input CHARGE="
            line_buffer.append(line)
            line = input_file_handle.readline()
            assert "SCANS=" in line, "Error: wrong input SCANS="
            line_buffer.append(line)
            #change the title with seq
            line = input_file_handle.readline()
            assert "RTINSECONDS=" in line, "Error: wrong input RTINSECONDS="
            line_buffer.append(line)
            line = input_file_handle.readline()
            assert "SEQ=" in line, "Error: wrong input SEQ="
            seq = line.split('=')[1]
            line_buffer.append(line)
            line_buffer[1] = f"TITLE={seq}"
            #add the sequence int the precursor dictionary
            # write the metadata in int ouput file
            for l in line_buffer:
              output_handle.write(l)
            while line and not ("END IONS" in line):
              line = input_file_handle.readline()
              output_handle.write(line)
            output_handle.write("\n")
            line = input_file_handle.readline()
            line = input_file_handle.readline()

def dict_to_pkl(fname, dict):
  with open (fname, mode = 'wb') as handle:
    pickle.dump( dict, handle)

def read_from_pkl(fname):
  d ={}
  with open(fname, mode='rb') as rhandle:
    d = pickle.load(rhandle)
  return d

def feature_list_for_annotation(input_spectrum_handle, spectrum_location_dict, spectrum_rtinseconds_dict, feature_lists):
    """TODO(nh2tran): docstring."""

    # ~ print("".join(["="] * 80)) # section-separating line
    # ~ print("WorkerIO: _parse_spectrum()")

    #read each feature from feature file
    for scan_list, rt_mean, precursor_mass, raw_sequence in feature_lists:
      neighbor_size = 5
      ### select best neighbors from the scan_list by their distance to rt_mean
      # probably move this selection to get_location(), run once rather than repeating
      neighbor_count = len(scan_list)
      best_scan_index = None
      best_distance = float('inf')
      for scan_index, scan in enumerate(scan_list):
        try:
          distance = abs(spectrum_rtinseconds_dict[scan] - rt_mean)
        except:
          distance = abs(spectrum_rtinseconds_dict[scan + '\n'] - rt_mean)

        if distance < best_distance:
          best_distance = distance
          best_scan_index = scan_index
      neighbor_center = best_scan_index
      neighbor_left_count = neighbor_center
      neighbor_right_count = neighbor_count - neighbor_left_count - 1
      neighbor_size_half = neighbor_size // 2
      neighbor_left_count = min(neighbor_left_count, neighbor_size_half)
      neighbor_right_count = min(neighbor_right_count, neighbor_size_half)



      ### parse and add neighbor spectra
      scan_list_middle = []

      for index in range(neighbor_center - neighbor_left_count, neighbor_center + neighbor_right_count + 1):
        scan = scan_list[index]
        scan_list_middle.append(scan)
      concat_label = '_'.join(scan_list_middle)

    return

def open_input(input_spectrum_file, input_feature_file):
  print("".join(["="] * 80))  # section-separating line
  print("WorkerIO: open_input()")
  input_spectrum_handle = open(input_spectrum_file, 'r')
  input_feature_handle = open(input_feature_file, 'r')
  return input_spectrum_handle, input_feature_handle

def get_location(input_spectrum_file):
  spectrum_location_file = input_spectrum_file + '.locations.pkl'
  if os.path.exists(spectrum_location_file):
    print("WorkerIO: read cached spectrum locations")
    with open(spectrum_location_file, 'rb') as fr:
      data = pickle.load(fr)
      spectrum_location_dict, spectrum_rtinseconds_dict, spectrum_count = data
  else:
    print("WorkerIO: build spectrum location from scratch")
    spectrum_location_dict = {}
    spectrum_rtinseconds_dict = {}
    line = True
    while line:
      current_location = input_spectrum_handle.tell()
      line = input_spectrum_handle.readline()
      if "BEGIN IONS" in line:
        spectrum_location = current_location
      elif "SCANS=" in line:
        scan = re.split('=|\r\n', line)[1]
        spectrum_location_dict[scan] = spectrum_location
      elif "RTINSECONDS=" in line:
        rtinseconds = float(re.split('=|\r\n', line)[1])
        spectrum_rtinseconds_dict[scan] = rtinseconds
    spectrum_location_dict = spectrum_location_dict
    spectrum_rtinseconds_dict = spectrum_rtinseconds_dict
    spectrum_count = len(spectrum_location_dict)
    with open(spectrum_location_file, 'wb') as fw:
      pickle.dump((spectrum_location_dict, spectrum_rtinseconds_dict, spectrum_count), fw)
  return spectrum_location_dict, spectrum_rtinseconds_dict

def extract_precursor_5_middile_info_by_spec_group_id(feature_f, spectrum_rtinseconds_dict):
  'return dictionary of {scan:(middle_scan_list, mz, z, rt_mean, seq)}'
  feature_df = pd.read_csv(feature_f, usecols=[0, 1, 2, 3, 4, 5])

  # print(feature_df)
  d = defaultdict(list)
  # with open(spec_file, mode= 'rb') as spec_f_h:
  num_scans = 0
  scan_l = []
  scan_set = set()
  test_counter = 0
  for _, row in feature_df.iterrows():
    test_counter+=1
    spec_group_id, mz, z, rt_mean, seq, scans = row
    mass = mz * z - mass_H * z
    if mass > MZ_MAX:
        print('Skipping the precursor {} with mass {} higher than maz max{}'.format(seq, mass, MZ_MAX))
        continue
    # find the spec_group_id in the dictionary
    scan_list = scans.split(';')
    num_scans += len(scan_list)
    neighbor_size = 5
    ### select best neighbors from the scan_list by their distance to rt_mean
    # probably move this selection to get_location(), run once rather than repeating
    neighbor_count = len(scan_list)
    best_scan_index = None
    best_distance = float('inf')
    for scan_index, scan in enumerate(scan_list):
      try:
        distance = abs(spectrum_rtinseconds_dict[scan] - rt_mean)
      except:
        distance = abs(spectrum_rtinseconds_dict[scan + '\n'] - rt_mean)
      #print(distance)
      if distance < best_distance:
        best_distance = distance
        best_scan_index = scan_index
      neighbor_center = best_scan_index
      neighbor_left_count = neighbor_center
      neighbor_right_count = neighbor_count - neighbor_left_count - 1
      neighbor_size_half = neighbor_size // 2
      neighbor_left_count = min(neighbor_left_count, neighbor_size_half)
      neighbor_right_count = min(neighbor_right_count, neighbor_size_half)

    ### parse and add neighbor spectra
    scan_list_middle = []

    for index in range(neighbor_center - neighbor_left_count, neighbor_center + neighbor_right_count + 1):
      scan = scan_list[index]
      scan_list_middle.append(scan)
    d[spec_group_id].append((scan_list_middle, mz, z, rt_mean, seq))
    #print(seq)
    #if seq == 'C(+57.02)SPHLVLSALTSDNHGATYAFSGTHYWR': print(seq)
  #dict_to_pkl('precursor_middle_scans.pkl', d)
  #print(test_counter, len(d))
  #print(d)
  return d

def extract_precursor_5_highest_intensity_info_by_spec_group_id(feature_f, input_spectrum_handle, spectrum_location_dict ):
  'return dictionary of {scan:(middle_scan_list, mz, z, rt_mean, seq)}'
  print("extract_precursor_5_highest_intensity_info_by_spec_group_id is calling ......")
  feature_df = pd.read_csv(feature_f, usecols=[0, 1, 2, 3, 4, 5])

  d = defaultdict(list)
  test_counter = 0
  for _, row in feature_df.iterrows():
    test_counter+=1
    spec_group_id, mz, z, rt_mean, seq, scans = row
    # find the spec_group_id in the dictionary
    scan_list = scans.split(';')
    scan_int = {}
    for scan in scan_list:
        try:
            spectrum_location = spectrum_location_dict[scan]
        except:
            spectrum_location = spectrum_location_dict[scan+'\n']
        input_spectrum_handle.seek(spectrum_location)
        # parse header lines
        line = input_spectrum_handle.readline()
        assert "BEGIN IONS" in line, "Error: wrong input BEGIN IONS"
        line = input_spectrum_handle.readline()
        assert "TITLE=" in line, "Error: wrong input TITLE="
        line = input_spectrum_handle.readline()
        assert "PEPMASS=" in line, "Error: wrong input PEPMASS="
        line = input_spectrum_handle.readline()
        assert "CHARGE=" in line, "Error: wrong input CHARGE="
        line = input_spectrum_handle.readline()
        assert "SCANS=" in line, "Error: wrong input SCANS="
        line = input_spectrum_handle.readline()
        assert "RTINSECONDS=" in line, "Error: wrong input RTINSECONDS="
        # parse fragment ions
        mz_list, intensity_list = _parse_spectrum_ion(input_spectrum_handle)
        normal_it = normalize_int(intensity_list)
        scan_int[scan] = np.sum(normal_it)
    top_five_int = dict(sorted(scan_int.items(), key=lambda item: item[1], reverse=True)[:5])
    ### parse and add neighbor spectra
    scan_list = list(top_five_int.keys())


    d[spec_group_id].append((scan_list, mz, z, rt_mean, seq))
  dict_to_pkl('precursor_middle_scans.pkl', d)
  return d

def extract_all_features_by_scan(feature_f, spectrum_rtinseconds_dict):
    'return dictionary of {scan:(middle_scan_list, mz, z, rt_mean, seq)}'
    feature_df = pd.read_csv(feature_f, usecols=[0, 1, 2, 3, 4, 5])
    d = defaultdict(list)
    num_scans = 0
    for _, row in feature_df.iterrows():
        spec_group_id, mz, z, rt_mean, seq, scans = row
    	# find the spec_group_id in the dictionary
        scan_list = scans.split(';')
        for scan in scan_list:
           d[scan].append((mz, z, rt_mean, seq))
        #d[spec_group_id].append((scan_list, mz, z, rt_mean, seq))
    return d

    
def extract_all_features_by_spec_group_id(feature_f, spectrum_rtinseconds_dict):
    'return dictionary of {scan:(middle_scan_list, mz, z, rt_mean, seq)}'
    feature_df = pd.read_csv(feature_f, usecols=[0, 1, 2, 3, 4, 5])
    d = defaultdict(list)
    num_scans = 0
    scan_l = []
    for _, row in feature_df.iterrows():
        spec_group_id, mz, z, rt_mean, seq, scans = row
    	# find the spec_group_id in the dictionary
        scan_list = scans.split(';')
        num_scans += len(scan_list)
        d[spec_group_id].append((scan_list, mz, z, rt_mean, seq))
    
    return d

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate MGF files for use in PepNet or Transformer_DIA.")
    parser.add_argument("--model_name", required=True, help="Name of the model being used for annotation (e.g., pepnet).")
    parser.add_argument("--spectrum_file", required=True, help="Path to the input spectrum .mgf file.")
    parser.add_argument("--feature_file", required=True, help="Path to the feature file (.csv).")
    parser.add_argument("--selection", required=True, choices=["all", "five_rt", "five_int"],
                        help="Annotation mode: 'all', 'five_rt' (top 5 closest to mean RT), or 'five_int' (top 5 intensity).")

    args = parser.parse_args()

    # Validate input files
    if not os.path.isfile(args.spectrum_file) or not args.spectrum_file.endswith(".mgf"):
        raise ValueError("Spectrum file must exist and have a .mgf extension.")
    if not os.path.isfile(args.feature_file) or not args.feature_file.endswith(".csv"):
        raise ValueError("Feature file must exist and have a .csv extension.")

    print("Annotating MGF file using model:", args.model_name)
    print("Spectrum file:", args.spectrum_file)
    print("Feature file:", args.feature_file)
    print("Annotation mode:", args.selection)

    input_spectrum_handle, input_feature_handle = open_input(args.spectrum_file, args.feature_file)
    spectrum_location_dict, spectrum_rtinseconds_dict = get_location(args.spectrum_file)

    precursor_info = {}

    if args.selection == 'all':
        all_precursor_info = extract_all_features_by_spec_group_id(args.feature_file, spectrum_rtinseconds_dict)
        annotate_mgf_one_spectrum_one_precursor(args.spectrum_file, spectrum_location_dict, all_precursor_info, args.model_name, True)

    elif args.selection == 'five_int':
        precursor_info = extract_precursor_5_highest_intensity_info_by_spec_group_id(args.feature_file, input_spectrum_handle, spectrum_location_dict)

    elif args.selection == 'five_rt':
        precursor_info = extract_precursor_5_middile_info_by_spec_group_id(args.feature_file, spectrum_rtinseconds_dict)

    else:
        sys.exit("Error: Invalid annotation mode.")


    out_name = args.feature_file.split('.csv')[0] + '_' + args.model_name + '_' + args.selection + '.mgf'
    annotate_mgf_one_spectrum_multiple_precursor(args.spectrum_file, spectrum_location_dict, precursor_info, args.model_name, out_name, args.selection)

