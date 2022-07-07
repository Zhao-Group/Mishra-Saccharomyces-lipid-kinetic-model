# Python3 code to integrate complete lipidomics pipeline for data analysis
# of yeast samples

import os, sys
import pandas as pd
import numpy as np
from statistics import mean
from read_input_file import read_input_file, read_input_file_2
from scan_averaging import scan_averaging_3
from offline_cal import offline_cal
from deisotope import deisotoping, deisotoping_13C_isotopes
from sterol_module import sterols_prm_product_ints

path = sys.argv[1]


# Reading the input file with m/z scan window information
pos_cycles, neg_cycles, pos_start_mz, pos_end_mz, neg_start_mz, neg_end_mz, pos_cal_peaks, neg_cal_peaks, database_filename, ppm_tol = read_input_file_2('input_settings_lipidomics.txt')


# Reading the lipid database
lipids_df = pd.ExcelFile('/'.join(['Database', database_filename])).parse()
lipid_dict = lipids_df.set_index('Lipid name').T.to_dict(orient='list')
sorted_dict = {k:v for k,v in sorted(lipid_dict.items(), key=lambda item: item[1][4])}      # (Sorting key, val) pairs by m/z value of monoisotopic M0 peak

# Sterol precursor:product m/z dictionary
precursor_mz = [411.4332, 414.3736]
product_mz = [376.3961, 379.3365]
sterol_dict = {int(key):val for key,val in zip(precursor_mz, product_mz)}       # Dictionary of precursor:product pairs for sterols
sterol_names = {411: 'IS Cholesterol-d7', 414: 'Ergosterol'}

# Handling Pos files:
files = os.listdir('/'.join([path, 'Pos']))
files.sort()
files = [file for file in files if os.path.splitext(file)[-1]=='.mzML']

for i in range(pos_cycles):
    detected_lipids = pd.DataFrame()
    print("Pos window: {0}".format(i+1))
    for file in files:
        avg_peaks = scan_averaging_3('/'.join([path, 'Pos', file]), pos_start_mz[i], pos_end_mz[i])                  # Scan averaging
        corr_peaks = offline_cal(avg_peaks, pos_cal_peaks[int(pos_start_mz[i])])                    # Offline calibration
        identified_peaks = deisotoping(corr_peaks, sorted_dict, ppm_tol, '+')                       # Deisotoping & Lipid identification
        detected_series = pd.Series({k:v for k,v in identified_peaks}, name=os.path.splitext(file)[0])
        detected_lipids = detected_lipids.append(detected_series)
    detected_lipids = detected_lipids.T.fillna(0).sort_index().astype('float')
    detected_lipids.to_excel('/'.join([path, 'Pos_' + str(int(pos_start_mz[i])) + '.xlsx']))

# Sterol handling
detected_lipids = pd.DataFrame()
for file in files:
    prod_ints = sterols_prm_product_ints('/'.join([path, 'Pos', file]), sterol_dict, ppm_tol)
    detected_series = pd.Series({sterol_names[key]:mean(prod_ints[key]) for key in prod_ints}, name=os.path.splitext(file)[0])
    detected_lipids = detected_lipids.append(detected_series)
detected_lipids = detected_lipids.T.fillna(0).sort_index().astype('float')
detected_lipids.to_excel('/'.join([path, 'Pos_400' + '.xlsx']))


# Handling Neg files:
files = os.listdir('/'.join([path, 'Neg']))
files.sort()
files = [file for file in files if os.path.splitext(file)[-1]=='.mzML']

for i in range(neg_cycles):
    detected_lipids = pd.DataFrame()
    print("Neg window: {0}".format(i+1))
    for file in files:
        print(file)
        avg_peaks = scan_averaging_3('/'.join([path, 'Neg', file]), neg_start_mz[i], neg_end_mz[i])                  # Scan averaging
        corr_peaks = offline_cal(avg_peaks, neg_cal_peaks[int(neg_start_mz[i])])                    # Offline calibration
        identified_peaks = deisotoping(corr_peaks, sorted_dict, ppm_tol, '-')                       # Deisotoping & Lipid identification
        detected_series = pd.Series({k:v for k,v in identified_peaks}, name=os.path.splitext(file)[0])
        detected_lipids = detected_lipids.append(detected_series)
    detected_lipids = detected_lipids.T.fillna(0).sort_index().astype('float')
    detected_lipids.to_excel('/'.join([path, 'Neg_' + str(int(neg_start_mz[i])) + '.xlsx']))
