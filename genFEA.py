## This function reads the vtu and vtp files from the image configuration
## and also the results and output the displacement field and vtk structures
## for the calculation and generation of the new FEA files

import numpy as np
import glob
from FluidVolume import FluidVolume
import IFEA_io as Iio
from vtkmodules.util import numpy_support as vtknp

def genFEA(fname, reference_mesh_path, f_ind, file_sub, p_info):
    '''
    This function reads the vtu and vtp files from the imaged configuration
    and also the results and output the displacement field and vtk structures
    for the calculation and generation of the new FEA files

    Args:
        - fname: str, the path to the files
        - f_ind: int, the index of the current simulation
        - file_sub: list, 
        - p_info: list, the list of landmark points
    '''

    if f_ind > 0:
        k = '%02d' % (f_ind - 1)
        fname_result_pre = "mesh_" + k + "1/"
    k = '%02d' % (f_ind)
    fname_result = "mesh_" + k + "result/"

    # Get file names of all mesh-complete files
    file_endo_lv = reference_mesh_path + "mesh-surfaces/endo_lv.vtp"
    file_endo_rv = reference_mesh_path+ "mesh-surfaces/endo_rv.vtp"
    file_top = reference_mesh_path + "mesh-surfaces/base.vtp"
    file_epi = reference_mesh_path + "mesh-surfaces/epi.vtp"
    file_epi_apex = reference_mesh_path + "mesh-surfaces/epi_apex.vtp"
    file_epi_mid = reference_mesh_path + "mesh-surfaces/epi_mid.vtp"
    file_vl = reference_mesh_path + "mesh-complete.mesh.vtu"
    folder_file = glob.glob(fname + fname_result + "*.vtu")
    folder_file.sort()

    # Get landmark points
    # These should be GlobalNodeIDs - 1 (0-indexed)
    lv_p_gindex = p_info[0]
    rv_p_gindex = p_info[1]
    epi_p_gindex = p_info[2]

    ## lv
    lv_output, lv_gnode = Iio.read_surf(file_endo_lv)
    lv_NoC = lv_output.GetNumberOfCells()
    conn_lv = []
    for i in range(lv_NoC):
        cell = lv_output.GetCell(i)
        npts = cell.GetNumberOfPoints()
        connt = [cell.GetPointId(j) for j in range(npts)]
        conn_lv.append(connt)
    conn_lv = np.asarray(conn_lv)

    ## rv
    rv_output, rv_gnode = Iio.read_surf(file_endo_rv)
    rv_NoC = rv_output.GetNumberOfCells()
    conn_rv = []
    for i in range(rv_NoC):
        cell = rv_output.GetCell(i)
        npts = cell.GetNumberOfPoints()
        connt = [cell.GetPointId(j) for j in range(npts)]
        conn_rv.append(connt)
    conn_rv = np.asarray(conn_rv)

    ## top
    top_output, top_gnode = Iio.read_surf(file_top)

    lv_edge_gnode = np.intersect1d(top_gnode, lv_gnode)
    rv_edge_gnode = np.intersect1d(top_gnode, rv_gnode)

    ## epi
    epi_output, epi_gnode = Iio.read_surf(file_epi)

    ## epi_apex
    epi_apex_output, epi_apex_gnode = Iio.read_surf(file_epi_apex)

    ## epi_mid
    epi_mid_output, epi_mid_gnode = Iio.read_surf(file_epi_mid)

    ## origin
    ref_output = Iio.read_vol(file_vl)
    point0d_vtk_array = ref_output.GetPoints().GetData()
    point0d = vtknp.vtk_to_numpy(point0d_vtk_array)

    # LV volume and landmark point positions
    lv_vlm_dat = FluidVolume(point0d, lv_edge_gnode, lv_gnode, conn_lv)
    lv_p_dat = np.zeros((len(lv_p_gindex), 3))
    for i in range(len(lv_p_gindex)):
        lv_p_dat[i, :] = point0d[lv_p_gindex[i], :]
    # RV volume and landmark point positions
    rv_vlm_dat = FluidVolume(point0d, rv_edge_gnode, rv_gnode, conn_rv)
    rv_p_dat = np.zeros((len(rv_p_gindex), 3))
    for i in range(len(rv_p_gindex)):
        rv_p_dat[i, :] = point0d[rv_p_gindex[i], :]
    # epi landmark point positions
    epi_p_dat = np.zeros((len(epi_p_gindex), 3))
    for i in range(len(epi_p_gindex)):
        epi_p_dat[i, :] = point0d[epi_p_gindex[i], :]

    ## ref+cur
    lv_vlm_cur = np.zeros(len(file_sub))
    rv_vlm_cur = np.zeros(len(file_sub))
    lv_p_cur = [None] * len(file_sub)
    rv_p_cur = [None] * len(file_sub)
    epi_p_cur = [None] * len(file_sub)
    Ind_result: int
    for Ind_result in range(len(file_sub)):
        file_rt = folder_file[round(len(folder_file) * file_sub[Ind_result]) - 1]
        cur_output = Iio.read_vol(file_rt)
        ## ref
        if Ind_result == 0:
            point0_vtk_array = cur_output.GetPoints().GetData()
            point0 = vtknp.vtk_to_numpy(point0_vtk_array)
            # LV
            lv_vlm_ref = FluidVolume(point0, lv_edge_gnode, lv_gnode, conn_lv)
            lv_p_ref = np.zeros((len(lv_p_gindex), 3))
            for i in range(len(lv_p_gindex)):
                lv_p_ref[i, :] = point0[lv_p_gindex[i], :]
            # RV
            rv_vlm_ref = FluidVolume(point0, rv_edge_gnode, rv_gnode, conn_rv)
            rv_p_ref = np.zeros((len(rv_p_gindex), 3))
            for i in range(len(rv_p_gindex)):
                rv_p_ref[i, :] = point0[rv_p_gindex[i], :]
            # epi
            epi_p_ref = np.zeros((len(epi_p_gindex), 3))
            for i in range(len(epi_p_gindex)):
                epi_p_ref[i, :] = point0[epi_p_gindex[i], :]
        ## cur
        disp_vtk_array = cur_output.GetPointData().GetArray('Displacement')
        disp = vtknp.vtk_to_numpy(disp_vtk_array)
        point1 = point0 + disp
        if Ind_result == 0:
            point1R = point1
        # LV
        lv_vlm_cur[Ind_result] = FluidVolume(point1, lv_edge_gnode, lv_gnode, conn_lv)
        lv_p_cur[Ind_result] = np.zeros((len(lv_p_gindex), 3))
        for i in range(len(lv_p_gindex)):
            lv_p_cur[Ind_result][i, :] = point1[lv_p_gindex[i], :]
        # RV
        rv_vlm_cur[Ind_result] = FluidVolume(point1, rv_edge_gnode, rv_gnode, conn_rv)
        rv_p_cur[Ind_result] = np.zeros((len(rv_p_gindex), 3))
        for i in range(len(rv_p_gindex)):
            rv_p_cur[Ind_result][i, :] = point1[rv_p_gindex[i], :]
        # epi
        epi_p_cur[Ind_result] = np.zeros((len(epi_p_gindex), 3))
        for i in range(len(epi_p_gindex)):
            epi_p_cur[Ind_result][i, :] = point1[epi_p_gindex[i], :]
    # print(lv_vlm_dat, lv_vlm_ref, lv_vlm_cur)
    # print(rv_vlm_dat, rv_vlm_ref, rv_vlm_cur)

    figure_vl_flag = False
    if figure_vl_flag == True:
        lv_vlm_cur_all = np.zeros(len(folder_file))
        rv_vlm_cur_all = np.zeros(len(folder_file))
        lv_p_cur_all = [None] * len(folder_file)
        rv_p_cur_all = [None] * len(folder_file)
        epi_p_cur_all = [None] * len(folder_file)
        for Ind_result in range(len(folder_file)):
            file_rt = folder_file[Ind_result]
            cur_output = Iio.read_vol(file_rt)
            disp_vtk_array = cur_output.GetPointData().GetArray('Displacement')
            disp = vtknp.vtk_to_numpy(disp_vtk_array)
            point1 = point0 + disp
            # LV
            lv_vlm_cur_all[Ind_result] = FluidVolume(point1, lv_edge_gnode, lv_gnode, conn_lv)
            lv_p_cur_all[Ind_result] = np.zeros((len(lv_p_gindex), 3))
            for i in range(len(lv_p_gindex)):
                lv_p_cur_all[Ind_result][i, :] = point1[lv_p_gindex[i], :]
            # RV
            rv_vlm_cur_all[Ind_result] = FluidVolume(point1, rv_edge_gnode, rv_gnode, conn_rv)
            rv_p_cur_all[Ind_result] = np.zeros((len(rv_p_gindex), 3))
            for i in range(len(rv_p_gindex)):
                rv_p_cur_all[Ind_result][i, :] = point1[rv_p_gindex[i], :]
            # epi
            epi_p_cur_all[Ind_result] = np.zeros((len(epi_p_gindex), 3))
            for i in range(len(epi_p_gindex)):
                epi_p_cur_all[Ind_result][i, :] = point1[epi_p_gindex[i], :]
        print(lv_vlm_dat, lv_vlm_ref, lv_vlm_cur_all)
        print(rv_vlm_dat, rv_vlm_ref, rv_vlm_cur_all)

    vals = {
        'point0d': point0d,
        'point0': point0,
        'point1R': point1R,
        'lv_output': lv_output,
        'lv_gnode': lv_gnode,
        'rv_output': rv_output,
        'rv_gnode': rv_gnode,
        'top_output': top_output,
        'top_gnode': top_gnode,
        'epi_output': epi_output,
        'epi_gnode': epi_gnode,
        'epi_apex_output': epi_apex_output,
        'epi_apex_gnode': epi_apex_gnode,
        'epi_mid_output': epi_mid_output,
        'epi_mid_gnode': epi_mid_gnode,
        'ref_output': ref_output,
        'cur_output': cur_output,
        'lv_vlm_cur': lv_vlm_cur,           # LV volume in the current configuration (cur)
        'rv_vlm_cur': rv_vlm_cur,           # RV volume in the current configuration (cur)
        'lv_vlm_dat': lv_vlm_dat,           # LV volume in the imaged/relaxed configuration (dat)
        'lv_vlm_ref': lv_vlm_ref,           # LV volume in the reference configuration (ref)
        'rv_vlm_dat': rv_vlm_dat,           # RV volume in the imaged/relaxed configuration (dat)
        'rv_vlm_ref': rv_vlm_ref,           # RV volume in the reference configuration (ref)
        'lv_edge_gnode': lv_edge_gnode,
        'rv_edge_gnode': rv_edge_gnode,
        'lv_p_dat': lv_p_dat,               # Position of LV landmark points in the imaged/relaxed configuration (dat)
        'rv_p_dat': rv_p_dat,               # Position of RV landmark points in the imaged/relaxed configuration (dat)
        'epi_p_dat': epi_p_dat,             # Position of epicardial landmark points in the imaged/relaxed configuration (dat)
        'lv_p_ref': lv_p_ref,               # Position of LV landmark points in the reference configuration (ref)
        'rv_p_ref': rv_p_ref,               # Position of RV landmark points in the reference configuration (ref)
        'epi_p_ref': epi_p_ref,             # Position of epicardial landmark points in the reference configuration (ref)
        'lv_p_cur': lv_p_cur,               # Position of LV landmark points in the current configuration (cur)
        'rv_p_cur': rv_p_cur,               # Position of RV landmark points in the current configuration (cur)
        'epi_p_cur': epi_p_cur              # Position of epicardial landmark points in the current configuration (cur)

    }
    return vals