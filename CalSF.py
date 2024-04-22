## This function provides a method to do optimization of the material parameters
## based on the current configuration and image configuration

import numpy as np

def CalSF(fname, f_ind, vals, imageData):
    '''
    This function calculates the stress free configuration based on the current
    configuration and image configuration

    Args:
        - fname: str, the path to the files
        - f_ind: int, the index of the current simulation
        - vals: dict, the dictionary of processed data (volumes, landmarks positions, etc.)
                from the most recent FEA simulation
        - imageData: dict, the dictionary of processed data from the image configuration (UNUSED)
    '''
    ####################### Calculation  #######################################################

    point0d = vals['point0d']
    # point0 = vals['point0']
    point1R = vals['point1R']

    if f_ind > 0:
        k = '%02d' % (f_ind - 1)
        fname_result_pre = "mesh_" + k + "1/"
    k = '%02d' % (f_ind)
    fname_result = "mesh_" + k + "result/"

    r_disp = np.zeros(len(point0d))
    r_pos = point1R - point0d
    for i in range(len(r_pos)):
        r_disp[i] = np.linalg.norm(r_pos[i, :])
    r_disp_max = np.max(r_disp)
    np.save(fname + fname_result + "r_disp.npy", r_pos)
    print("Max disp diff: " + str(r_disp_max))

    # beta = 1
    # if "fname_result_pre" in globals():
    #     r_pos_1 = np.load(fname + fname_result_pre + "r_disp.npy")
    #     beta = -beta * np.tensordot(r_pos_1, (r_pos - r_pos_1), axes=2) / np.tensordot((r_pos - r_pos_1),
    #                                                                                    (r_pos - r_pos_1), axes=2)
    # print(beta)
    beta = 1

    lv_vlm_cur = vals['lv_vlm_cur']
    rv_vlm_cur = vals['rv_vlm_cur']
    print('lv_vlm_dat: ', vals['lv_vlm_dat'], 'lv_vlm_ref: ', vals['lv_vlm_ref'], 'lv_vlm_cur: ', lv_vlm_cur)
    print('rv_vlm_dat: ', vals['rv_vlm_dat'], 'rv_vlm_ref: ', vals['rv_vlm_ref'], 'rv_vlm_cur: ', rv_vlm_cur)

    ####################### Finish Calculation  #######################################################
    cal_vals = {
        'beta': beta,
        'r_pos': r_pos,
        'r_disp_max': r_disp_max,
    }

    return cal_vals
