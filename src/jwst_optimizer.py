#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Practical for course 'Space Instrumentation',
TU Delft, The Netherlands
2022
By L. Welzel, Z. Burr
"""

import numpy as np
from scipy.optimize import minimize
from os import getcwd
import matlab.engine
from static_plotters import displC, plotCAmpl

def plot_psf(correction_matrix, show_MATLAB=False):
    out = eng.JWST_Simulation(matlab.double(correction_matrix.tolist()), show_MATLAB, nargout=12)

    psf = np.array(out[7]._data).reshape(out[7].size, order='F').T
    psf_diff = np.array(out[8]._data).reshape(out[8].size, order='F').T
    displC(np.log10(psf), title=r"Corrected PSF (log$_{10}$-scale)", trim=251)
    displC(np.log10(psf_diff), title=r"Diffraction limited PSF (log$_{10}$-scale)", trim=251)

    displC(np.log10(np.abs(psf - psf_diff)), title=r"PSF difference (log$_{10}$-scale)", trim=251)



def objective(theta):
    theta = theta.reshape(matrix_shape)
    # (opd_rms, strehl, spotsize_rms, mtf) = eng.JWST_sim_runtime(matlab.double(theta.tolist()), nargout=4)
    out = eng.JWST_sim_runtime(matlab.double(theta.tolist()), nargout=4)
    # out is (opd_rms, strehl, spotsize_rms, mtf)
    out = np.array(out) / expected_values * minimization_rescaling
    out[1] = 1 / out[1]
    out = np.mean(np.square(out))
    print(out)
    if np.allclose(np.clip(out, a_min=1., a_max=None), 1., atol=1e-01):
        return 1.
    return out


def optimize_jwst_mirror_segments(correction_mat):
    # plot_psf(correction_mat, show_MATLAB=False)

    print(f"\t Starting optimization...\n"
          f"\t\t Optimizing {correction_mat.shape[1]} parameters of {correction_mat.shape[0]} reflector segments...")

    result = minimize(objective, x0=correction_mat,
                      method="SLSQP",
                      options={"disp": True,
                               "maxiter": 50,
                               "ftol": 1.e-6})
    print(f"\tConverged educed chi-squared: {result.fun}")
    print(f"\tFunction evaluations: {result.nfev}\n"
          f"\t\t The comparison of models with a reduced chi-squared between\n"
          f"\t\t{1 - np.sqrt(2/result.nfev) * 3:.3f} and "
          f"\t\t{1 + np.sqrt(2/result.nfev) * 3:.3f} (3 sigma) is uncertain.")

    converged = result.x.reshape(matrix_shape)
    plot_psf(converged, show_MATLAB=True)
    return


if __name__ == "__main__":
    global eng
    global matrix_shape
    global minimization_rescaling
    global expected_values
    print("Starting up MATLAB engine...")
    eng = matlab.engine.start_matlab()
    eng.addpath(eng.genpath(getcwd()), nargout=0)
    print("MATLAB engine running.")

    matrix_shape = (18, 4)

    ### EXPECTED VALUES FOR CHI_SQURED OPTIMIZATION
    softness_scale = 1
    minimization_rescaling = softness_scale * np.array([1., 1., 1., 1.])
    expected_opd = 0.01  # 0.0746 -> in practice diff limited
    expected_strehl = 1.
    expected_spotsize_rms = 0.003
    expected_mtf = 1.

    expected_values = np.array([expected_opd, expected_strehl, expected_spotsize_rms, expected_mtf])

    correction_mat = np.array([
        [0.414170806677252, 0.300243140031183, 0.439976576044261, 0.0152119184255706, -0.128208842154579,
         0.478345378762405, 0.0754076647493227, 0.651252356266327, 0.101404608428473, 0.0391465366416755,
         0.548089404530392, 0.266178356301017, -0.360151894856671, 0.385371897751932, 0.300917231394459,
         -0.546668246684149, -0.238745289802985, -0.585858177909011],

        [0.139276147276064, -0.390138294249314, -0.411334159967717, 0.338255587537226, 0.350357337374621,
         -0.244629702055557, -0.188140054852467, -0.208542872352273, 0.252227970049943, 0.0246373453963109,
         0.411647424007853, 0.210408685278169, 0.472651076977497, -0.182166946273771, 0.0847186192633207,
         -0.390257631406097, -0.230116336295599, -0.0477925462370182],

        [-0.214891914341358, -0.0239202817325437, 0.326573979042765, 0.0709910754624060, -0.438971070807491,
         -0.309013559302602, 0.327732173448263, -0.292396965620019, 0.176871093438419, -0.181895273849737,
         -0.352344222848263, -0.330232933973511, 0.448108735396022, -0.366189014643874, 0.171462889478031,
         -0.0574700377971163, -0.106588493632423, 0.0846413033551107],

        [-4.77158933834332, -5.90414896183012, 6.46915662700510, 12.3803995282438, -7.83571389646577, 7.02863171336911,
         -7.60505047433170, 12.3089743294982, 7.55653924410722, -1.89433225241735, 5.72159670128435,
         -0.0411604369112084, 7.72475667997612, -3.58727666297385, -10.6689141377548, 2.27478638187123,
         10.2546957682035, -7.65585159583474]
    ])

    rng = np.random.default_rng()
    noise = rng.normal(0, 0.0001, correction_mat.shape) * np.mean(np.abs(correction_mat), axis=0)

    correction_mat += noise

    correction_mat = correction_mat.T

    optimize_jwst_mirror_segments(correction_mat)

    input("Press any key to continue...\n"
          "\t This will close the MATLAB figures.")

    eng.quit()
