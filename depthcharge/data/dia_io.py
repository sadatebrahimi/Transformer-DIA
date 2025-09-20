import numpy as np


# Special vocabulary symbols - we always put them at the start.
_PAD = "_PAD"
_GO = "_GO"
_EOS = "_EOS"
_START_VOCAB = [_PAD, _GO, _EOS]

PAD_ID = 0
GO_ID = 1
EOS_ID = 2
vocab_reverse = ['A',
                 'R',
                 'N',
                 'N(Deamidation)',
                 'D',
                 #~ 'C',
                 'C(Carbamidomethylation)',
                 'E',
                 'Q',
                 'Q(Deamidation)',
                 'G',
                 'H',
                 'I',
                 'L',
                 'K',
                 'M',
                 'M(Oxidation)',
                 'F',
                 'P',
                 'S',
                 'T',
                 'W',
                 'Y',
                 'V',
                ]

vocab_reverse = _START_VOCAB + vocab_reverse

vocab = dict([(x, y) for (y, x) in enumerate(vocab_reverse)])

vocab_size = len(vocab_reverse)
mass_H = 1.0078
mass_H2O = 18.0106
mass_NH3 = 17.0265
mass_N_terminus = 1.0078
mass_C_terminus = 17.0027
mass_CO = 27.9949


mass_AA = {'_PAD': 0.0,
           '_GO': mass_N_terminus-mass_H,
           '_EOS': mass_C_terminus+mass_H,
           'A': 71.03711, # 0
           'R': 156.10111, # 1
           'N': 114.04293, # 2
           'N(Deamidation)': 115.02695,
           'D': 115.02694, # 3
           #~ 'C(Carbamidomethylation)': 103.00919, # 4
           'C(Carbamidomethylation)': 160.03065, # C(+57.02)
           #~ 'C(Carbamidomethylation)': 161.01919, # C(+58.01) # orbi
           'E': 129.04259, # 5
           'Q': 128.05858, # 6
           'Q(Deamidation)': 129.0426,
           'G': 57.02146, # 7
           'H': 137.05891, # 8
           'I': 113.08406, # 9
           'L': 113.08406, # 10
           'K': 128.09496, # 11
           'M': 131.04049, # 12
           'M(Oxidation)': 147.0354,
           'F': 147.06841, # 13
           'P': 97.05276, # 14
           'S': 87.03203, # 15
           'T': 101.04768, # 16
           'W': 186.07931, # 17
           'Y': 163.06333, # 18
           'V': 99.06841, # 19
          }

mass_ID = [mass_AA[vocab_reverse[x]] for x in range(vocab_size)]
mass_ID_np = np.array(mass_ID, dtype=np.float32)

mass_AA_min = mass_AA["G"] # 57.02146
num_ion = 8
class WorkerDIA(object):
    def __init__(self):
        self.MZ_MAX = 3000
        self.SPECTRUM_RESOLUTION = 50 # bins for 1.0 Da = precision 0.02 Da
        self.MZ_SIZE = int(self.MZ_MAX * self.SPECTRUM_RESOLUTION)
        self.neighbor_size = 5
        self.mass_H = 1.0078
        self.dia_window = 10
        self._buckets =[12, 22, 32]
        self.vocab = vocab
        self.GO_ID =GO_ID
        self.EOS_ID = EOS_ID
        self.PAD_ID = PAD_ID
        self.vocab_size = vocab_size
        self.WINDOW_SIZE = self.dia_window
        self.SPECTRUM_RESOLUTION = self.SPECTRUM_RESOLUTION
        self.mass_ID_np = mass_ID_np
        self.mass_H2O = mass_H2O
        self.mass_NH3 = mass_NH3
        self.num_ion = num_ion
        self.mass_C_terminus = mass_C_terminus
        self.mass_N_terminus = mass_N_terminus
        self.mass_ID = mass_ID

    def process_spectrum(self, spectrum_mz_list, spectrum_intensity_list, peptide_mass):


        # neutral mass, location, assuming ion charge z=1
        charge = 1.0
        spectrum_mz = np.array(spectrum_mz_list, dtype=np.float32)
        neutral_mass = spectrum_mz - charge * self.mass_H
        neutral_mass_location = np.rint(neutral_mass * self.SPECTRUM_RESOLUTION).astype(
            np.int32)
     
        neutral_mass_location_view = neutral_mass_location

        # intensity
        spectrum_intensity = np.array(spectrum_intensity_list, dtype=np.float32)
        spectrum_intensity_max = np.max(spectrum_intensity)
        norm_intensity = spectrum_intensity
        norm_intensity_view = norm_intensity

        # fill spectrum holders
        spectrum_holder = np.zeros(shape=(1, self.MZ_SIZE), dtype=np.float32)

        spectrum_holder_view = spectrum_holder

        for index in range(neutral_mass_location.size):
            spectrum_holder_view[0, neutral_mass_location_view[index]] = max(
                spectrum_holder_view[0, neutral_mass_location_view[index]],  
                norm_intensity_view[index])  
        spectrum_original_forward = np.copy(spectrum_holder)
        spectrum_original_backward = np.copy(spectrum_holder)


        # peptide_mass
        spectrum_original_forward[0, int(round(
            peptide_mass * self.SPECTRUM_RESOLUTION))] = spectrum_intensity_max
        spectrum_original_backward[0, int(round(
            peptide_mass * self.SPECTRUM_RESOLUTION))] = spectrum_intensity_max
        # N-terminal, b-ion, peptide_mass_C
        # append N-terminal
        
        # append peptide_mass_C
        mass_C = self.mass_C_terminus + self.mass_H
        peptide_mass_C = peptide_mass - mass_C
        spectrum_original_forward[0, int(round(
            peptide_mass_C * self.SPECTRUM_RESOLUTION))] = spectrum_intensity_max
        # C-terminal, y-ion, peptide_mass_N
        
        # append peptide_mass_N
        mass_N = self.mass_N_terminus - self.mass_H
        peptide_mass_N = peptide_mass - mass_N
        spectrum_original_backward[0, int(round(
            peptide_mass_N * self.SPECTRUM_RESOLUTION))] = spectrum_intensity_max

        return spectrum_original_forward, spectrum_original_backward

    def get_location(self, peptide_mass, prefix_mass, direction):
        if direction == 0:
            candidate_b_mass = prefix_mass + self.mass_ID_np
            candidate_y_mass = peptide_mass - candidate_b_mass
        elif direction == 1:
            candidate_y_mass = prefix_mass + self.mass_ID_np
            candidate_b_mass = peptide_mass - candidate_y_mass

        # b-ions
        candidate_b_H2O = candidate_b_mass - self.mass_H2O
        candidate_b_NH3 = candidate_b_mass - self.mass_NH3
        candidate_b_plus2_charge1 = ((candidate_b_mass + 2 * self.mass_H) / 2
                                     - self.mass_H)

        # y-ions
        candidate_y_H2O = candidate_y_mass - self.mass_H2O
        candidate_y_NH3 = candidate_y_mass - self.mass_NH3
        candidate_y_plus2_charge1 = ((candidate_y_mass + 2 * self.mass_H) / 2
                                     - self.mass_H)
        # ion_8
        b_ions = [candidate_b_mass,
                  candidate_b_H2O,
                  candidate_b_NH3,
                  candidate_b_plus2_charge1]
        y_ions = [candidate_y_mass,
                  candidate_y_H2O,
                  candidate_y_NH3,
                  candidate_y_plus2_charge1]
        ion_mass_list = b_ions + y_ions
        ion_mass = np.array(ion_mass_list, dtype=np.float32)

        # ion locations
        location_sub50 = np.rint(ion_mass * self.SPECTRUM_RESOLUTION).astype(np.int32)
        location_sub50 -= (self.WINDOW_SIZE // 2)
        location_plus50 = location_sub50 + self.WINDOW_SIZE
        ion_id_rows, aa_id_cols = np.nonzero(np.logical_and(
            location_sub50 >= 0,
            location_plus50 <= self.MZ_SIZE))
        return ion_id_rows, aa_id_cols, location_sub50, location_plus50

    
    def copy_values(self, candidate_intensity_view, spectrum_view, location_sub,i1, i2):
    
        i1_start = self.neighbor_size * i1
        for neighbor in range(self.neighbor_size):
            for j in range(self.WINDOW_SIZE):
                try:
                    candidate_intensity_view[i2, i1_start + neighbor, j] = spectrum_view[neighbor, location_sub[i1, i2] + j]
                except:
                    continue
                    
    def get_candidate_intensity(self, spectrum_original, peptide_mass, prefix_mass, direction):
        ion_id_rows, aa_id_cols, location_sub50, location_plus50 = self.get_location(peptide_mass, prefix_mass, direction)
        # candidate_intensity
        candidate_intensity = np.zeros(shape=(self.vocab_size,
                                          self.neighbor_size * self.num_ion,
                                          self.WINDOW_SIZE),
                                   dtype=np.float32)
   
        location_sub50_view = location_sub50
 
        location_plus50_view = location_plus50
        candidate_intensity_view = candidate_intensity
        row = ion_id_rows.astype(np.int32)
        col = aa_id_cols.astype(np.int32)
        
        for index in range(ion_id_rows.size):
            if col[index] < 3:
                continue
            self.copy_values(candidate_intensity_view, spectrum_original, location_sub50_view, row[index], col[index])
        max_intensity = np.max(candidate_intensity)
        if max_intensity > 1.0:
            candidate_intensity /= max_intensity
        return candidate_intensity

    def _parse_spectrum(self, precursor_mass, mz_lists, intensity_lists, neighbor_right_count, neighbor_size_half):
        
        MZ_SIZE = int(self.MZ_MAX * self.SPECTRUM_RESOLUTION)
        spectrum_original_forward_list = []
        spectrum_original_backward_list = []
        for mz_list, intensity_list in zip(mz_lists, intensity_lists):

            spectrum_original_forward, spectrum_original_backward = self.process_spectrum(mz_list,
                                                                 intensity_list,
                                                                 precursor_mass)
            spectrum_original_forward_list.append(spectrum_original_forward)
            spectrum_original_backward_list.append(spectrum_original_backward)

        if neighbor_right_count < neighbor_size_half:
            for x in range(neighbor_size_half - neighbor_right_count):
                spectrum_original_forward_list.append(np.zeros(
                    shape=(1, self.MZ_SIZE),
                    dtype=np.float32))
                spectrum_original_backward_list.append(np.zeros(
                    shape=(1, self.MZ_SIZE),
                    dtype=np.float32))

        spectrum_original_forward = np.vstack(spectrum_original_forward_list)
        spectrum_original_backward = np.vstack(spectrum_original_backward_list)
       
        return spectrum_original_forward, spectrum_original_backward

    def _process_raw_seq(self, raw_seq):
        peptide = []
        raw_len = len(raw_seq)
        pep_len, index = 0, 0
        while index < raw_len:
            if raw_seq[index] == '+':
                if peptide[-1] == 'C' and raw_seq[index:index + 7] == '+57.021':
                    peptide[-1] = 'C(Carbamidomethylation)'
                    index += 7
            else:
                peptide.append(raw_seq[index])
                index += 1
                pep_len += 1
        return peptide, pep_len

    def calculate_ms2(self, precursor_mz, precursor_charge, scan_list_middle, mz_list, int_list, neighbor_right_count, neighbor_size_half, seq):
        candidate_intensity_list_forward = []
        candidate_intensity_list_backward = []

        precursor_mass = precursor_mz * precursor_charge - self.mass_H * precursor_charge

        if precursor_mass > self.MZ_MAX:
            return None
        spectrum_original_forward, spectrum_original_backward = self._parse_spectrum(precursor_mass, mz_list, int_list, neighbor_right_count, neighbor_size_half)
        peptide, peptide_len = self._process_raw_seq(seq)
        for bucket_id, target_size in enumerate(self._buckets):
            if peptide_len + 2 <= target_size:  # +2 to include GO and EOS
                break
        decoder_size = self._buckets[bucket_id]
        peptide_ids = [self.vocab[x] for x in peptide]
        pad_size = decoder_size - (len(peptide_ids) + 2)
        # forward

        peptide_ids_forward = peptide_ids[:]
        peptide_ids_forward.insert(0, self.GO_ID)
        peptide_ids_forward.append(self.EOS_ID)
        peptide_ids_forward += [self.PAD_ID] * pad_size
        
        peptide_ids_backward = peptide_ids[::-1]
        peptide_ids_backward.insert(0, self.EOS_ID)
        peptide_ids_backward.append(self.GO_ID)
        peptide_ids_backward += [self.PAD_ID] * pad_size

        prefix_mass = 0.0
        suffix_mass = 0.0
        for index in range(self._buckets[-1]):
            if index < decoder_size:
                prefix_mass += self.mass_ID[peptide_ids_forward[index]]
                candidate_intensity_forward = self.get_candidate_intensity(
                    spectrum_original_forward,
                    precursor_mass,
                    prefix_mass,
                    0)
                cand_int_shape_forward = candidate_intensity_forward.shape

                suffix_mass += self.mass_ID[peptide_ids_backward[index]]
                candidate_intensity_backward = self.get_candidate_intensity(
                    spectrum_original_backward,
                    precursor_mass,
                    suffix_mass,
                    1)
                cand_int_shape_backward = candidate_intensity_backward.shape
            else:
                candidate_intensity_forward = np.zeros(cand_int_shape_forward)
                candidate_intensity_backward = np.zeros(cand_int_shape_backward)

            candidate_intensity_list_forward.append(candidate_intensity_forward)
            candidate_intensity_list_backward.append(candidate_intensity_backward)
        return candidate_intensity_list_forward, candidate_intensity_list_backward
