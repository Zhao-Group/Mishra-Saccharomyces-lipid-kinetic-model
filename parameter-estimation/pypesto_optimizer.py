import os, sys
from os.path import join
import pandas as pd
import numpy as np
import random
import math

import amici
import petab
import libsbml

import pypesto
import pypesto.optimize as optimize
import pypesto.petab

# os.environ['BLAS_LIBS'] = '-lopenblas'


petab_folder = sys.argv[1]
results_folder = sys.argv[2]

model_name = 'lipid_model'

path_to_folder = '/home/smishr10/amici_pyPESTO'
yaml_config = join(path_to_folder, 'petab_files', petab_folder, 'petab_problem.yaml')

if not os.path.isdir(join(path_to_folder, 'results', results_folder)):
    os.mkdir(join(path_to_folder, 'results', results_folder))

if not os.path.isdir(join(path_to_folder, 'results', results_folder, 'history')):
    os.mkdir(join(path_to_folder, 'results', results_folder, 'history'))


petab_problem = petab.Problem.from_yaml(yaml_config)
importer = pypesto.petab.PetabImporter(petab_problem, model_name=model_name)


print('\n\n---------------Starting model creation----------------\n\n')

model = importer.create_model()
solver = importer.create_solver(model)


objective = importer.create_objective(model, solver, guess_steadystate=False)    

objective.amici_solver.setSensitivityMethod(amici.SensitivityMethod.adjoint)
objective.amici_solver.setSensitivityOrder(amici.SensitivityOrder.second)
objective.amici_solver.setRelativeTolerance(1e-3)
objective.amici_model.requireSensitivitiesForAllParameters()


n_starts = 1000

problem = importer.create_problem(objective)


options = {'maxiter': 50, 'disp': False}
optimizer = optimize.ScipyOptimizer(options=options)



engine = pypesto.engine.MultiProcessEngine()

hist_opt = pypesto.objective.HistoryOptions(trace_record=True, trace_save_iter=10,
        storage_file=join(path_to_folder, 'results', results_folder, 'history', 'history_{id}.csv'))


result = optimize.minimize(problem=problem, optimizer=optimizer, n_starts=n_starts, engine=engine, history_options=hist_opt,
        filename=join(path_to_folder, 'results', results_folder, 'result_file.hdf5'))

print('Results of fvals from optimization')
print(result.optimize_result.get_for_key("fval"))

df = result.optimize_result.as_dataframe()

df.to_excel(join(path_to_folder, 'results', results_folder, 'test_optimizer_result.xlsx'))
