##### cpjohn expansion microscope

### Last Update
# Ver. 1.1
# Date: 2020/10/22
# Content: Modify layout, add dtype function.

### Before
# Ver. 1.0
# R. Date : 2020/10/09

import numpy as np
from scipy.interpolate import RegularGridInterpolator

def tag_to_posi(tag):
    temp_posi = list(map(int, tag.split('-')))
    return temp_posi

def posi_to_tag(posi):
    temp_tag = '-'.join(list(map(str, posi)))
    return temp_tag

def ReLU(x):
    y = 0
    if x > 0: y = x
    return y

def AutoPC(mx_c, mx_t, posi_c, posi_t): # Proto in 10 T-07
    # Function ReLU
    def ReLU(x):
        y = 0
        if x > 0: y = x
        return y
    # diff_a : Distance from posi_c anchor to posi_t anchor
    diff_a = [t-c for t, c in zip(posi_t, posi_c)]
    # Get anchor distance along corresponding axis
    dist_z, dist_y, dist_x = diff_a
    # Get side length of matrices
    mxc_z, mxc_y, mxc_x = np.shape(mx_c)
    mxt_z, mxt_y, mxt_x = np.shape(mx_t)
    # Trimming overlap area
    mx_c_overlap = mx_c[ReLU( dist_z) : mxc_z          - ReLU(mxc_z - mxt_z - dist_z),
                        ReLU( dist_y) : mxc_y          - ReLU(mxc_y - mxt_y - dist_y),
                        ReLU( dist_x) : mxc_x          - ReLU(mxc_x - mxt_x - dist_x)]
    mx_t_overlap = mx_t[ReLU(-dist_z) : mxc_z - dist_z - ReLU(mxc_z - mxt_z - dist_z),
                        ReLU(-dist_y) : mxc_y - dist_y - ReLU(mxc_y - mxt_y - dist_y),
                        ReLU(-dist_x) : mxc_x - dist_x - ReLU(mxc_x - mxt_x - dist_x)]
    # Get matrix total item amount
    overlap_box_size = np.prod(np.shape(mx_c_overlap))
    mx_c_overlap_1d = np.reshape(mx_c_overlap, overlap_box_size)
    mx_t_overlap_1d = np.reshape(mx_t_overlap, overlap_box_size)
    PC = np.corrcoef(mx_c_overlap_1d, mx_t_overlap_1d)[0][1]
    return PC

def resampleRGI3d(input_mx, resize_to, dtype='float64'):
    # Ref: 10 3-A-7
    # Ver. 2.1
    # input_mx : numpy array, the original target matrix
    # resize_to: list or tuple with 3 int inside
    a, b, c = np.shape(input_mx)
    p, q, r = resize_to
    z_grid = np.linspace(0, p - 1, a)
    y_grid = np.linspace(0, q - 1, b)
    x_grid = np.linspace(0, r - 1, c)
    RGI = RegularGridInterpolator((z_grid, y_grid, x_grid), input_mx)
    z_grid_t2 = np.arange(p)
    y_grid_t2 = np.arange(q)
    x_grid_t2 = np.arange(r)
    meshgrid_para = np.meshgrid(z_grid_t2, y_grid_t2, x_grid_t2)
    RGI_mesh_mx = RGI((meshgrid_para[0], meshgrid_para[1], meshgrid_para[2]))
    RGI_mx = np.transpose(RGI_mesh_mx, axes=[1, 0, 2]).astype(dtype)
    return RGI_mx