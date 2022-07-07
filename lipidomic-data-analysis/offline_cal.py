# Python3 code to perform offline calibration on scan-averaged peak lists

import numpy as np

def error(a,b):
    return abs((a - b) * 1E6/a)

def return_mz(raw_peaks, cal_peaks, ppm_tol):
    cal_peaks_in_data = []
    peaks_in_data = []
    for mz in cal_peaks:
        for mz_raw in raw_peaks:
            if(error(mz, mz_raw) <= ppm_tol):
                cal_peaks_in_data.append(mz)
                peaks_in_data.append(mz_raw)
                continue
    return cal_peaks_in_data, peaks_in_data


def offline_cal(peaks, cal_peaks):
    cal_peaks_in_data, peaks_in_data = return_mz(peaks[:,0], cal_peaks, 10)         # Allowing an error tolerance of 10 ppm between real and exact m/z
    # print(cal_peaks_in_data)
    # print(peaks_in_data)
    pol = np.polyfit(peaks_in_data, cal_peaks_in_data, 1)
    offline_cal_peaks = np.polyval(pol, peaks[:,0])
    peaks[:,0] = offline_cal_peaks

    return peaks


if __name__ == '__main__':
    print()
