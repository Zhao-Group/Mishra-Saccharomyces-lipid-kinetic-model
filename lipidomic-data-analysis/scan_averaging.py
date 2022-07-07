# Python3 function to perform scan averaging of MS spectra from lipidomics experiments
# stored in .mzML files

import sys
import pymzml
import numpy as np
from math import sqrt, fabs
import pandas as pd


def bin_max_mz(bin_start, m_R, R_orbitrap):
    R_real = R_orbitrap*sqrt(m_R/bin_start)
    return bin_start*(1.0 + 1.0/R_real)

def mz_I_avg_1st(bin):
    I_avg = np.mean([i[1] for i in bin])
    mz_avg = sum(i[0] * i[1] for i in bin)/(len(bin) * I_avg)
    all_mzs = [i[0] for i in bin]
    all_Is = [i[1] for i in bin]
    p_n = [mz_avg, I_avg, all_mzs, all_Is]
    return p_n

def mz_I_avg_later(bin):
    I_avg = np.mean([i[1] for i in bin])
    mz_avg = sum(i[0] * i[1] for i in bin)/(len(bin) * I_avg)
    all_mzs = [i for p in bin for i in p[2]]
    all_Is = [i for p in bin for i in p[3]]
    p_n = [mz_avg, I_avg, all_mzs, all_Is]
    return p_n

def mz_I_avg(bin, iter):
    I_avg = np.mean([i[1] for i in bin])
    mz_avg = sum(i[0] * i[1] for i in bin)/(len(bin) * I_avg)
    if iter == 0:
            all_mzs = [i[0] for i in bin]
            all_Is = [i[1] for i in bin]
    else:
            all_mzs = [i for p in bin for i in p[2]]
            all_Is = [i for p in bin for i in p[3]]
    p_n = [mz_avg, I_avg, all_mzs, all_Is]
    return p_n

def mz_I_avg_final(peak):
    I_avg = np.mean([i for i in peak[3]])
    mz_avg = sum(i[0] * i[1] for i in zip(peak[2], peak[3]))/(len(peak[2]) * I_avg)
    return mz_avg, I_avg

# clustered_peaks = []
# clustered_intensities = []
# count = 0
# # val = all_peaks[0]
# cluster = []
# cluster_int = []
# int_cov = []
# while (count < len(sorted_peaks)-1):
#     cluster.append(sorted_peaks[count])
#     cluster_int.append(sorted_intensities[count])
#     if(error(mean(cluster),sorted_peaks[count+1]) > 10.0):
#         if(len(cluster) > 1):
#             int_cov.append(stdev(cluster_int)*100/mean(cluster_int))
#             print('Mean m/z: {:.4f}\tError m/z interval: {:.1f}\tMean Int: {:.2E}\tInt COV: {:.2f}%'.format(mean(cluster), error(mean(cluster), mean(cluster)-stdev(cluster)), mean(cluster_int), int_cov[-1]))
#             clustered_peaks.append(mean(cluster))
#             clustered_intensities.append(mean(cluster_int))
#         cluster = []
#         cluster_int = []
#         val = all_peaks[count+1]
#     count += 1


def scan_averaging(filename, start_mz):
    msrun = pymzml.run.Reader(filename)

    m_R = 200
    R_orbitrap = 140000

    # peaks_concatenate = np.concatenate([msrun[i].peaks('centroided') for i in range(1,msrun.get_spectrum_count()+1) if ((msrun[i].ms_level == 1) and (abs(min(msrun[i].mz) - start_mz) <= 50))])      # Slower, buggier implementation

    peaks_concatenate = np.concatenate([spec.peaks('centroided') for spec in msrun if spec.ms_level == 1 and fabs(min(spec.mz) - start_mz) <= 50])
    # timeit.timeit(lambda: np.concatenate([msrun[i].peaks('centroided') for i in range(1,msrun.get_spectrum_count()+1) if ((msrun[i].ms_level == 1) and (abs(min(msrun[i].mz) - 500) <= 50))]), number=1)

    # Sort peaks_concatenate (which is of shape (N,2)) by the column of index 0
    peaks_sort = peaks_concatenate[peaks_concatenate[:,0].argsort()]
    peaks_list = [[i[0], i[1], []] for i in peaks_sort]

    for iter in range(0,3):
        count = 0
        Bin = [peaks_list[0]]
        peaks_new = []
        bin_max = bin_max_mz(Bin[0][0], m_R, R_orbitrap)
        while (count < len(peaks_list)):
            if(peaks_list[count][0] < bin_max):
                Bin.append(peaks_list[count])
            else:
                if (iter == 0):
                    peaks_new.append(mz_I_avg_1st(Bin))
                else:
                    peaks_new.append(mz_I_avg_later(Bin))
                Bin = [peaks_list[count]]
                bin_max = bin_max_mz(Bin[0][0], m_R, R_orbitrap)
            count += 1
        # if (iter == 0):
        #     peaks_new.append(mz_I_avg_1st(Bin))
        # else:
        #     peaks_new.append(mz_I_avg_later(Bin))
        peaks_list = peaks_new

    Spec_final = []
    for p in peaks_list:
        Spec_final.append(mz_I_avg(p))

    Spec_final = np.asarray(Spec_final)
    return Spec_final


def scan_averaging_2(filename):
    msrun = pymzml.run.Reader(filename)

    m_R = 200
    R_orbitrap = 140000

    # peaks_concatenate = np.concatenate([msrun[i].peaks('centroided') for i in range(1,msrun.get_spectrum_count()+1) if ((msrun[i].ms_level == 1) and (abs(min(msrun[i].mz) - start_mz) <= 50))])      # Slower, buggier implementation

    peaks_concatenate = np.concatenate([spec.peaks('centroided') for spec in msrun if spec.ms_level == 1])
    # timeit.timeit(lambda: np.concatenate([msrun[i].peaks('centroided') for i in range(1,msrun.get_spectrum_count()+1) if ((msrun[i].ms_level == 1) and (abs(min(msrun[i].mz) - 500) <= 50))]), number=1)

    # Sort peaks_concatenate (which is of shape (N,2)) by the column of index 0
    peaks_sort = peaks_concatenate[peaks_concatenate[:,0].argsort()]
    peaks_list = [[i[0], i[1], []] for i in peaks_sort]

    for iter in range(0,3):
        count = 0
        Bin = [peaks_list[0]]
        peaks_new = []
        bin_max = bin_max_mz(Bin[0][0], m_R, R_orbitrap)
        while (count < len(peaks_list)):
            if(peaks_list[count][0] < bin_max):
                Bin.append(peaks_list[count])
            else:
                if (iter == 0):
                    peaks_new.append(mz_I_avg_1st(Bin))
                else:
                    peaks_new.append(mz_I_avg_later(Bin))
                Bin = [peaks_list[count]]
                bin_max = bin_max_mz(Bin[0][0], m_R, R_orbitrap)
            count += 1
        # if (iter == 0):
        #     peaks_new.append(mz_I_avg_1st(Bin))
        # else:
        #     peaks_new.append(mz_I_avg_later(Bin))
        peaks_list = peaks_new

    Spec_final = []
    for p in peaks_list:
        Spec_final.append(mz_I_avg(p))

    Spec_final = np.asarray(Spec_final)
    return Spec_final


def scan_averaging_3(filename, start_mz, end_mz):
    msrun = pymzml.run.Reader(filename)

    m_R = 200
    R_orbitrap = 140000

    peaks_concatenate = np.concatenate([spec.peaks('centroided') for spec in msrun if spec.ms_level == 1 and fabs(min(spec.mz) - start_mz) <= 15 and fabs(max(spec.mz) - end_mz) <= 15])

    # timeit.timeit(lambda: np.concatenate([msrun[i].peaks('centroided') for i in range(1,msrun.get_spectrum_count()+1) if ((msrun[i].ms_level == 1) and (abs(min(msrun[i].mz) - 500) <= 50))]), number=1)

    # Sort peaks_concatenate (which is of shape (N,2)) by the column of index 0
    peaks_sort = peaks_concatenate[peaks_concatenate[:,0].argsort()]
    peaks_list = [[i[0], i[1], [], []] for i in peaks_sort]

    for iter in range(0,3):
        Bin = [peaks_list[0]]
        bin_max = bin_max_mz(Bin[0][0], m_R, R_orbitrap)
        peaks_new = []
        for peak in peaks_list[1:]:
            if peak[0] < bin_max:
                Bin.append(peak)
            else:
                peaks_new.append(mz_I_avg(Bin, iter))
                Bin = [peak]
                bin_max = bin_max_mz(Bin[0][0], m_R, R_orbitrap)
        peaks_new.append(mz_I_avg(Bin, iter))
        peaks_list = peaks_new

    Spec_final = []
    for p in peaks_list:
        Spec_final.append(mz_I_avg_final(p))

    Spec_final = np.asarray(Spec_final)
    return Spec_final

if __name__ == '__main__':
    filename = sys.argv[1]
    start_mz = sys.argv[2]
    end_mz = sys.argv[3]
    Spec_final = scan_averaging_3(filename, float(start_mz), float(end_mz))
    # print(Spec_final)
    df = pd.DataFrame(Spec_final, columns=['m/z', 'Intensity'])
    df.to_excel('temp.xlsx', index=False)
