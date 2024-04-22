'''This code is designed to find both the stress-free configuration
and the optimized material parameters for the passive process
from the diastasis to the end-diastolic states
HO model is used as the constitutive model and GA is used for
the IFEA process
by Lei Shi at Vedula Lab in Columbia University 2022
'''

import numpy as np
import os
import shutil
from runFEA import runFEA
from getIMG import getIMG
from CalSF import CalSF
from Calculation import Calculation
from genFEA import genFEA
from genTar import genTar
from scipy.optimize import least_squares
import GeAlgo as ga
import time



# ----------------------- Global parameters -------------------------------

# Path to mesh-complete folder (mesh in relaxed/imaged configuration)
reference_mesh_path = "/Users/aaronbrown/Library/CloudStorage/GoogleDrive-abrown97@stanford.edu/My Drive/Stanford/Marsden Lab/Papers/Patient-specific coupled BiV Paper/Sims/mesh_and_fibers/meshes/truncated_BiV-mesh-complete/"

# Command to run svFSI simulation
exec_svfsi = "mpiexec -np 4 /Users/aaronbrown/Documents/GitHub/svFSI_vvedula22/build/svFSI-build/bin/svFSI "


def Simulations(x):
    '''
    Given material parameters x, this function determines the unloaded/reference
    configuration, using Sellier's method. It then computes the objective function
    based on the target volumes and control point displacements.
    
    This function executes one iteration of the outer optimization loop.
    '''

    # Copy reference mesh to mesh_00.
    if os.path.isdir("mesh_00") == False:
        shutil.copytree(reference_mesh_path, "mesh_00")

    # Get material parameters
    m_a = 10 ** x[0]    # Convert from log10 to normal
    m_b = x[1]
    m_af = 10 ** x[2]
    m_bf = x[3]
    m_as = 10 ** x[4]
    m_bs = x[5]
    m_afs = 2160
    m_bfs = 11.25
    m_eta = 400 # 400
    m_E = 1.0e5

    # Store material parameters (optimized and fixed)
    mat_para_opt = np.array([m_a, m_b, m_af, m_bf, m_as, m_bs, m_afs, m_bfs])
    mat_para_fix = np.array([m_eta, m_E])

    # Node indices on reference mesh for selected landmarks
    # TODO: Modify these to be the correct indices. These should be GlobalNodeID-1
    lv_p_gindex = np.array([3299,2755,413])
    rv_p_gindex = np.array([657,2836,396])
    epi_p_gindex = np.array([]) # Not needed
    p_info = [lv_p_gindex, rv_p_gindex, epi_p_gindex]

    # Path to reference meshes warped by image data
    fname = os.getcwd() + "/"
    imagePath = "../truncated_BiV-mesh-complete_morph/"
    file_index = ["RR70", "RR80", "RR90", "RR0"]
    file_sub = np.array([0.25,0.5,0.75,1])

    r_disp_max = 1
    finalflag = False

    t = time.time()

    # Inner loop to find stress free configuration
    for f_ind in range(3):
        # f_ind = 3
        print("++++++++++++++++++++++++++++++++++++++")
        print("The " + str(f_ind) + " generation: ")

        runFEA(mat_para_opt, f_ind, mat_para_fix, exec_svfsi, finalflag)
        vals = genFEA(fname, reference_mesh_path, f_ind, np.array([1]), p_info)
        imageData = getIMG(imagePath, file_index, vals, p_info)
        cal_vals = CalSF(fname, f_ind, vals, imageData)
        genTar(fname, reference_mesh_path, f_ind, vals, cal_vals)
        r_disp_max = cal_vals['r_disp_max']
        if r_disp_max < 0.01:
            break


    print("++++++++++++++++++++++++++++++++++++++")
    print("The final " + str(f_ind) + " generation: ")
    finalflag = True
    runFEA(mat_para_opt, f_ind, mat_para_fix, exec_svfsi, finalflag)
    vals = genFEA(fname, reference_mesh_path, f_ind, file_sub, p_info)
    cal_vals = Calculation(fname, f_ind, vals, imageData)

    k = '%02d' % (f_ind + 1)
    shutil.rmtree("mesh_" + k, ignore_errors=True)
    print("Elapsed: " + str(time.time() - t) + " seconds")

    fit_err = cal_vals['fit_err']
    if np.isnan(fit_err):
        fit_err = 100

    return fit_err

# -------------------------- Main code ----------------------------------------#

# Navigate to the directory of the script
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Start timer
tinit = time.time()

# -------------------------- Genetic algorithm --------------------------------#
# Lower and upper bounds for parameters
# a, b, a_f, b_f, a_s, b_s NOTE: a, a_f, a_s are log10-ed to
# avoid order effects
lb = np.array([2,1,3,1,3,1])
ub = np.array([5,40,6,40,6,40])

algorithm_parameter = {'max_num_iteration': 30,
                        'population_size':4,
                        'mutation_probability':0.1,
                        'mutation_change_generation': 5,
                        'mutation_change_factor': 1.3}

r = ga.GeAlgo(Simulations, lb, ub, algorithm_parameter)
# -----------------------------------------------------------------------------#

# -------------------------------- Bayesian -----------------------------------#

# ----------------------------------------------------------------------------#

print("Total Elapsed: " + str(time.time() - tinit))

print()



