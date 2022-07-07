# Python code for sterol analysis from Pos_AP files

import sys
import pymzml
import numpy as np
import math
import pandas as pd
from offline_cal import error
from deisotope import find_peak


def sterols_prm_product_ints(filename, sterol_dict, ppm_tol):
    msrun = pymzml.run.Reader(filename)

    prod_ints = {int(key):[] for key in sterol_dict}

    for spec in msrun:
        if spec.ms_level == 2 and int(spec.selected_precursors[0]['mz']) in sterol_dict:
            prec = int(spec.selected_precursors[0]['mz'])
            prod = sterol_dict[prec]
            peaks = spec.peaks('centroided')
            idx = find_peak(peaks[:,0], prod)
            if error(prod, peaks[idx,0]) < ppm_tol:
                # print('Error: {}\tProduct ion: {}\tDetected peak: {}'.format(error(prod, peaks[idx,0]), prod, peaks[idx,0]))
                prod_ints[prec].append(peaks[idx,1])
    for key in prod_ints:
        if not prod_ints[key]: prod_ints[key] = [0]

    return prod_ints



if __name__ == '__main__':
    filename = sys.argv[1]
    ppm_tol = 5.0

    precursor_mz = [411.4332, 414.3736]
    product_mz = [376.3961, 379.3365]
    sterol_dict = {int(key):val for key,val in zip(precursor_mz, product_mz)}       # Dictionary of precursor:product pairs for sterols

    prod_ints = sterols_prm_product_ints(filename, sterol_dict, ppm_tol)
    # print(prod_ints)
