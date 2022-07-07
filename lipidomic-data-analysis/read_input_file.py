# Python3 function to read the input parameters of a Lipidomics pipeline
# from a .txt file

def read_input_file(input_file):
    with open(input_file, 'r') as f_in:
        d = {}
        for line in f_in:
            if(line[0] not in ('#','\n')):
                (key, val) = line.split(': ')
                d[key] = val

    pos_cycles = int(d['Pos Ion cycling windows'])
    neg_cycles = int(d['Neg Ion cycling windows'])

    pos_start_mz = [float(i) for i in d['Pos Starting m/z values of window'].split(', ')]
    pos_end_mz = [float(i) for i in d['Pos Ending m/z values of window'].split(', ')]

    neg_start_mz = [float(i) for i in d['Neg Starting m/z values of window'].split(', ')]
    neg_end_mz = [float(i) for i in d['Neg Ending m/z values of window'].split(', ')]

    pos_cal_peaks = {}
    pos_cal_peaks[220] = [float(i) for i in d['Pos Cal peaks 220'].split(', ')]
    pos_cal_peaks[500] = [float(i) for i in d['Pos Cal peaks 500'].split(', ')]

    neg_cal_peaks = {}
    neg_cal_peaks[200] = [float(i) for i in d['Neg Cal peaks 200'].split(', ')]
    neg_cal_peaks[505] = [float(i) for i in d['Neg Cal peaks 505'].split(', ')]

    database_filename = d['Lipids database'].rstrip()
    ppm_tol = float(d['Mass accuracy tolerance (ppm)'])

    return pos_cycles, neg_cycles, pos_start_mz, pos_end_mz, neg_start_mz, neg_end_mz, pos_cal_peaks, neg_cal_peaks, database_filename, ppm_tol


def read_input_file_2(input_file):
    with open(input_file, 'r') as f_in:
        d = {}
        for line in f_in:
            if(line[0] not in ('#','\n')):
                (key, val) = line.split(': ')
                d[key] = val

    pos_cycles = int(d['Pos Ion cycling windows'])
    neg_cycles = int(d['Neg Ion cycling windows'])

    pos_start_mz = [float(i) for i in d['Pos Starting m/z values of window'].split(', ')]
    pos_end_mz = [float(i) for i in d['Pos Ending m/z values of window'].split(', ')]

    neg_start_mz = [float(i) for i in d['Neg Starting m/z values of window'].split(', ')]
    neg_end_mz = [float(i) for i in d['Neg Ending m/z values of window'].split(', ')]

    pos_cal_peaks = {}
    pos_cal_peaks[420] = [float(i) for i in d['Pos Cal peaks 420'].split(', ')]
    pos_cal_peaks[500] = [float(i) for i in d['Pos Cal peaks 500'].split(', ')]
    pos_cal_peaks[734] = [float(i) for i in d['Pos Cal peaks 734'].split(', ')]

    neg_cal_peaks = {}
    neg_cal_peaks[350] = [float(i) for i in d['Neg Cal peaks 350'].split(', ')]
    neg_cal_peaks[545] = [float(i) for i in d['Neg Cal peaks 545'].split(', ')]
    neg_cal_peaks[640] = [float(i) for i in d['Neg Cal peaks 640'].split(', ')]
    neg_cal_peaks[790] = [float(i) for i in d['Neg Cal peaks 790'].split(', ')]

    database_filename = d['Lipids database'].rstrip()
    ppm_tol = float(d['Mass accuracy tolerance (ppm)'])

    return pos_cycles, neg_cycles, pos_start_mz, pos_end_mz, neg_start_mz, neg_end_mz, pos_cal_peaks, neg_cal_peaks, database_filename, ppm_tol
