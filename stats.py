import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import warnings
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm

asteroid = sys.argv[1]

d2r = np.pi / 180.
r2d = 180. / np.pi

def _p1(x):
    return np.exp(-90.56*np.square(np.tan(x/2)))*(1 - (0.986*np.sin(x))/(0.119 + 1.341*np.sin(x) - 0.754*np.square(np.sin(x)))) + (1-np.exp(-90.56*np.square(np.tan(x/2))))*np.exp(-3.332*np.power(np.tan(1/2 * x), 0.631))

def _p2(x):
    return np.exp(-90.56*np.square(np.tan(x/2)))*(1 - (0.238*np.sin(x))/(0.119 + 1.341*np.sin(x) - 0.754*np.square(np.sin(x)))) + (1-np.exp(-90.56*np.square(np.tan(x/2))))*np.exp(-1.862*np.power(np.tan(1/2 * x), 1.218))

def HG(alpha, H, G):
    return H - 2.5 * np.log10((1-G) * _p1(alpha) + G * _p2(alpha))

C_type = pd.read_csv(f'{asteroid}-C_type_models.txt', delimiter='\t', names=['ModelNo', 'PhaseAngle', 'ReducedMag'])
S_type = pd.read_csv(f'{asteroid}-S_type_models.txt', delimiter='\t', names=['ModelNo', 'PhaseAngle', 'ReducedMag'])

model_numbers_C = C_type.ModelNo.unique()
model_numbers_S = S_type.ModelNo.unique()

uncertainties = pd.read_csv(f'{asteroid}-geom.txt', delimiter='\t', header=None)[6].values

C_H_arr = np.zeros(shape=len(model_numbers_C))
C_G_arr = np.zeros(shape=len(model_numbers_C))
x_plot_range = np.linspace(0, C_type.PhaseAngle.max(), 100)
failed_models_C = []
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for i in tqdm(model_numbers_C):

        failed_to_fit = False
        model_data = C_type[C_type.ModelNo == i]

        try:
            (H, G), _ = curve_fit(HG, model_data.PhaseAngle.values, model_data.ReducedMag.values, sigma=uncertainties, absolute_sigma=True, bounds=([0, -1], [30, 1]))
        except RuntimeError:
            print(f"Model {i} could not fit. Skipping.")
            failed_to_fit = True
            failed_models_C.append(i)

        if not failed_to_fit:
            C_H_arr[i], C_G_arr[i] = H, G
        
C_H_std = np.std(C_H_arr)
C_G_std = np.std(C_G_arr)

S_H_arr = np.zeros(shape=len(model_numbers_S))
S_G_arr = np.zeros(shape=len(model_numbers_S))
x_plot_range = np.linspace(0, C_type.PhaseAngle.max(), 100)
failed_models_S = []
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for i in tqdm(model_numbers_S):
        
        failed_to_fit = False
        model_data = S_type[S_type.ModelNo == i]
        
        try:
            (H, G), _ = curve_fit(HG, model_data.PhaseAngle.values, model_data.ReducedMag.values, sigma=uncertainties, absolute_sigma=True, bounds=([0, -1], [30, 1]))
        except RuntimeError:
            print(f"Model {i} could not fit. Skipping.")
            failed_to_fit = True
            failed_models_S.append(i)

        if not failed_to_fit:
            S_H_arr[i], S_G_arr[i] = H, G
        
S_H_std = np.std(S_H_arr)
S_G_std = np.std(S_G_arr)

if failed_models_S or failed_models_C:
    failed_models_all = list(set(failed_models_S + failed_models_C))
    
    np.delete(C_H_arr, failed_models_all)
    np.delete(C_G_arr, failed_models_all)
    np.delete(S_H_arr, failed_models_all)
    np.delete(S_G_arr, failed_models_all)
    
H_uncs = [C_H_std, S_H_std]
G_uncs = [C_G_std, S_G_std]
weights = [0.476, 0.524]
H_aspect_uncertainty = np.average(H_uncs, weights=weights)
G_aspect_uncertainty = np.average(G_uncs, weights=weights)

print(f'H Aspect Uncertainty: {H_aspect_uncertainty:.2f} mag')
print(f'G Aspect Uncertainty: {G_aspect_uncertainty:.2f}')