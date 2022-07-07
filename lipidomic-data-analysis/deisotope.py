# Python3 code to run deisotoping on offline-calibrated m/z peaks; this code
# sequentially runs a Type II and a Type I isotopic correction

import numpy as np
import math
from offline_cal import error
from molmass import Formula
import re

# def temp1(array, value):                    # Slower peak finding function; no longer used
#     array = np.asarray(array)
#     idx = (np.abs(array - value)).argmin()
#     return array[idx]

def find_peak(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def isotopologues(string_formula):
    array = []
    stuff = re.findall('C[0-9]+', string_formula)[0]
    carbon_len = int(stuff.replace('C',''))
    remainder = string_formula.replace(stuff,'')

    for p in range(carbon_len+1):
        a = carbon_len - p
        b = p
        if a !=0 and b!= 0:
            iso_form = 'C' + str(a) + '[13C]' + str(b) + remainder
        elif b == 0:
            iso_form = 'C' + str(a) + remainder
        elif a == 0:
            iso_form = '[13C]' + str(b) + remainder
        f = Formula(iso_form)
        d = dict(f.spectrum())
        a = list(d.keys())
        b = []
        for p in a[:4]:
            b = b + d[p]
        array.append(b)
    return array




def deisotoping(peaks, lipids_db, ppm_tol, charge_sign):              # Sorted dictionary lipids_db is passed here
    identified_peaks = []
    for lipid in lipids_db:
        mz = lipids_db[lipid][4]
        idx = find_peak(peaks[:,0], mz)
        if error(mz, peaks[idx,0]) < ppm_tol and lipids_db[lipid][2][0] == charge_sign:
            type1_corr = lipids_db[lipid][5]
            M0_intensity = peaks[idx,1]
            identified_peaks.append([lipid, M0_intensity/type1_corr])       # Type I correction applied here
            # Apply Type II isotopic correction here
            isotopes = lipids_db[lipid][6:]
            for i in range(int(len(isotopes)/2)):
                mz_iso = isotopes[2*i]
                fraction_iso = isotopes[2*i + 1]
                idx_iso = find_peak(peaks[:,0], mz_iso)
                if error(mz_iso, peaks[idx_iso,0] < ppm_tol):
                    peaks[idx_iso,1] = max(peaks[idx_iso,1] - M0_intensity * fraction_iso/type1_corr, 0)
    return np.asarray(identified_peaks)


def identification(peaks, lipids_db, ppm_tol, charge_sign):              # Sorted dictionary lipids_db is passed here
    identified_peaks = []
    for lipid in lipids_db:
        mz = lipids_db[lipid][4]
        idx = find_peak(peaks[:,0], mz)
        if error(mz, peaks[idx,0]) < ppm_tol and lipids_db[lipid][2][0] == charge_sign:
            M0_intensity = peaks[idx,1]
            identified_peaks.append([lipid, M0_intensity])
    return np.asarray(identified_peaks)



def deisotoping_13C_isotopes(peaks, lipids_db, ppm_tol, charge_sign):              # Sorted dictionary lipids_db is passed here
    identified_peaks = []
    for lipid in lipids_db:
        isotopologue_array = isotopologues(lipids_db[lipid][3])
        for num, isotope in enumerate(isotopologue_array):
            mz = isotope[0]
            idx = find_peak(peaks[:,0], mz)
            if error(mz, peaks[idx,0]) < ppm_tol and lipids_db[lipid][2][0] == charge_sign:
                type1_corr = isotope[1]
                M0_intensity = peaks[idx,1]
                identified_peaks.append([lipid + ' : ' + str(num), M0_intensity/type1_corr])       # Type I correction applied here

                # Apply Type II isotopic correction here
                isotopes = isotope[2:]
                for i in range(int(len(isotopes)/2)):
                    mz_iso = isotopes[2*i]
                    fraction_iso = isotopes[2*i + 1]
                    idx_iso = find_peak(peaks[:,0], mz_iso)
                    if error(mz_iso, peaks[idx_iso,0] < ppm_tol):
                        peaks[idx_iso,1] = max(peaks[idx_iso,1] - M0_intensity * fraction_iso/type1_corr, 0)
    return np.asarray(identified_peaks)