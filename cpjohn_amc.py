##### cpjohn_amc.py
#####
##### Last Update: 20200101
#####
##### "amc" stands for AM file configureration.
#####
##### It's a manual module made up for Analysis of Drosophila brain neuron
##### 3-D matrix data. This py module file contents functions below: 
#####


import numpy as np
import matplotlib.pyplot as plt
import os
import re


def f_reader(file_path, readtype = 'r'):
    target_txt = open(file_path, '{}'.format(readtype))
    target_read = target_txt.read()
    target_txt.close()
    return target_read

def happy_time(start,stop):
    process_time = round(stop - start)
    ss = process_time % 60
    mm = process_time // 60 % 60
    hh = process_time // 3600
    duration = "Process time == {}s == {}H {}m {}s".format(process_time,hh,mm,ss)
    return duration

def eucdist(position1, position2):
    if len(position1) != len(position2): print('Fail: Dimension not match')
    if type(position1) != list: print('Fail: Input type should be list')
    if type(position2) != list: print('Fail: Input type should be list')
    dist2 = 0.
    for i in range(len(position1)):
        dist2 += (position1[i] - position2[i]) ** 2
    dist = np.sqrt(dist2)
    return dist

def angle_between(vector1, vector2):
    #if len(vector1) != len(vector2): print('Fail: Dimension not match')
    #if len(vector1) != 2 and len(vector1) != 3: print('Fail: vector1 length error')
    #if len(vector2) != 2 and len(vector2) != 3: print('Fail: vector2 length error')
    #if type(vector1) != list: print('Fail: Input type should be list')
    #if type(vector2) != list: print('Fail: Input type should be list')
    vec1_2 = 0.
    vec2_2 = 0.
    dot_product = 0.
    for i in range(len(vector1)):
        vec1_2 += vector1[i] ** 2
        vec2_2 += vector2[i] ** 2
        dot_product += vector1[i] * vector2[i]
    len_vec1 = np.sqrt(vec1_2)
    len_vec2 = np.sqrt(vec2_2)
    angle = np.arccos(dot_product / len_vec1 / len_vec2) / np.pi * 180
    return angle

def brtxtCP(brtxt_read, swc=False):
    spliter = '\n------------------------------\n'
    br_li= list(filter(None, brtxt_read.split('@ ')))
    total_cp_li = []
    for i in range(len(br_li)): #len(br_li)
        br_meta = list(filter(None, br_li[i].split(spliter)))
        if swc == False: # Original branches.txt route
            cp = [br_meta[1]] + br_meta[2].split('\n')
            cp_li = list(map(lambda x: list(map(int, x.split(' ')[:3][::-1])), cp))
        elif swc == True: # *.swc route
            cp = br_meta[1].split('\n')
            cp_li = list(map(lambda x: list(map(float, x.split(' ')[:3])), cp))
        else:
            print('Error: swc setting failed.')
            break
        total_cp_li.append(cp_li)
    return total_cp_li

def vector2(position1, position2):
    if len(position1) != len(position2): print('Fail: Dimension not match')
    if type(position1) != list: print('Fail: Input type should be list')
    if type(position2) != list: print('Fail: Input type should be list')
    vec2 = [y - x for (x, y) in zip(position1, position2)]
    return vec2

def mx_preview(array):
    #plt.figure(figsize=(8, 8))
    plt.imshow(array, aspect = 'equal', interpolation='nearest')
    plt.title('Preview of (x, y) = {}'.format(np.shape(array)),
              fontsize = 20, y=1.05)
    plt.xlabel('X', fontsize = 20)
    plt.ylabel('Y', fontsize = 20)
    plt.colorbar(fraction=0.04, pad=0.1)
    return plt.show()

def mx_preview2(array, min=0, max=2048, plotsize=6):
    plt.figure(figsize=(plotsize, plotsize))
    plt.imshow(array, aspect = 'equal', interpolation='nearest')
    plt.title('(x, y) = {}, Int limit = {}'.format(np.shape(array), [min, max]),
              fontsize = 16, y=1.05)
    plt.xlabel('X', fontsize = 16)
    plt.ylabel('Y', fontsize = 16)
    plt.clim(min,max)
    plt.colorbar(fraction=0.0451, pad=0.05)
    return plt.show()

def sth(train_history, train, validation): # sth stands for show_train_history
    plt.plot(train_history.history[train])
    plt.plot(train_history.history[validation])
    plt.title('Train History')
    plt.ylabel(train)
    plt.xlabel('Epoch')
    plt.legend(['train', 'validation'], loc = 'upper left')
    return plt.show()

def plot_images_labels_prediction(images, labels, prediction, idx, num = 30): # From BOOK
    fig = plt.gcf()
    fig.set_size_inches(12, 14)
    if num > 30: num = 30 
    for i in range(0, num):
        ax = plt.subplot(6, 5, 1+i)
        ax.imshow(images[idx],)
        title = "label=" + str(labels[idx])
        if len(prediction) > 0:
            title += ", predict=" + str(prediction[idx]) 
        ax.set_title(title, fontsize = 10) 
        ax.set_xticks([]); ax.set_yticks([])        
        idx += 1 
    plt.show()

def proj_xy_2d(path_of_am):
    target_am = open(path_of_am, 'rt', encoding = 'utf8')
    read_am = target_am.read()
    target_am.close()
    for item in read_am.split("\n"):
        if 'define Lattice' in item:
            key_lattice = item.strip()
        else:
            pass
    x_latti = int(key_lattice.split(' ')[-3])
    y_latti = int(key_lattice.split(' ')[-2])
    z_latti = int(key_lattice.split(' ')[-1])
    pre_li = list(map(lambda x: int(x), filter(None, read_am.split('@1\n')[2].split('\n'))))
    pre_ar1 = np.array(pre_li).reshape(z_latti, y_latti, x_latti)
    project_mx = np.amax(pre_ar1, axis = 0)
    return project_mx

####################################################################################################

"""
def soma_index(path_of_branches_txt):
    # Part of loading soma position (proto) in branch.txt
    brtxt_filename = str(path_of_branches_txt)
    target_bptxt = open(brtxt_filename, 'rt', encoding = 'utf8')
    read_bptxt = target_bptxt.read()[:1000]
    target_bptxt.close()
    trim_soma = re.findall('--*.\n[0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+\n--*.\n',read_bptxt)
    soma_index_proto = list(map(int, trim_soma[0].split('\n')[1].split(' ')[-5::-1]))
    return soma_index_proto
"""

"""
def ep_index(path_of_EPAM):
    EPAM_filename = str(path_of_EPAM)
    target_EPAM = open(EPAM_filename, 'rt', encoding = 'utf8')
    read_EPAM = target_EPAM.read()
    target_EPAM.close()
    trim_EP1 = read_EPAM.split('@')[-1].split('\n')
    if '' in trim_EP1: trim_EP1.remove('')
    if '1' in trim_EP1: trim_EP1.remove('1')
    trim_EP2 = list(map(lambda x: list(map(int,x.split(' '))), trim_EP1)) # 3D data
    return trim_EP2
"""

"""
def rdm_pick_slice_no(number, rdm_list):
    record_name = rdm_list[number]
    count = 0
    while rdm_list[number] == rdm_list[number-1]:
        count += 1
        number += -1
    return [record_name, count]
"""

"""
def soma_position(path_of_branches_txt):
    # Part of loading soma position (proto) in branch.txt
    brtxt_filename = str(path_of_branches_txt)
    target_bptxt = open(brtxt_filename, 'rt', encoding = 'utf8')
    read_bptxt = target_bptxt.read()[:1000]
    target_bptxt.close()
    trim_soma = re.findall('--*.\n[0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+\n--*.\n',read_bptxt)
    soma_posi_proto = list(map(int, trim_soma[0].split('\n')[1].split(' ')[-5::-1]))
    return soma_posi_proto
"""

"""
def cali_soma_position(path_of_branches_txt, path_of_sfAM):
    # Part of loading soma position (proto) in branch.txt
    brtxt_filename = str(path_of_branches_txt)
    target_bptxt = open(brtxt_filename, 'rt', encoding = 'utf8')
    read_bptxt = target_bptxt.read()[:1000]
    target_bptxt.close()
    trim_soma = re.findall('--*.\n[0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+\n--*.\n',read_bptxt)
    soma_posi_proto = list(map(int, trim_soma[0].split('\n')[1].split(' ')[-5::-1]))
    # Part of loading calibration axis in SF.am
    sfAM_filename = str(path_of_sfAM)
    target_sfAM = open(sfAM_filename, 'rt', encoding = 'utf8')
    read_sfAM = target_sfAM.read()[:1000]
    target_sfAM.close()
    grab_ca = re.findall('Min. [X,Y,Z] = [0-9]+', read_sfAM)
    calibrate_axis = list(map(lambda x : int(x.split(' = ')[-1]),grab_ca))
    # Part of merge
    soma_posi = [0,0,0]
    for m in range(3):
        soma_posi[m] = int(soma_posi_proto[m] - calibrate_axis[m] + 1)
    return soma_posi
"""

"""
def cali_EP_position(path_of_EPAM, path_of_SFAM):
    # Part of loading EP position
    EPAM_filename = str(path_of_EPAM)
    target_EPAM = open(EPAM_filename, 'rt', encoding = 'utf8')
    read_EPAM = target_EPAM.read()
    target_EPAM.close()
    trim_EP1 = read_EPAM.split('@')[-1].split('\n')
    if '' in trim_EP1: trim_EP1.remove('')
    if '1' in trim_EP1: trim_EP1.remove('1')
    trim_EP2 = list(map(lambda x: list(map(int,x.split(' '))), trim_EP1)) # 3D data
    # Part of loading calibration axis in SF.am
    sfAM_filename = path_of_SFAM
    target_sfAM = open(sfAM_filename, 'rt', encoding = 'utf8')
    read_sfAM = target_sfAM.read()[:1000]
    target_sfAM.close()
    grab_ca = re.findall('Min. [X,Y,Z] = [0-9]+', read_sfAM)
    calibrate_axis = list(map(lambda x : int(x.split(' = ')[-1]),grab_ca))
    # Part of merge
    for n in range(len(trim_EP2)):
        trim_EP2[n][0] = trim_EP2[n][0] - calibrate_axis[0] + 1
        trim_EP2[n][1] = trim_EP2[n][1] - calibrate_axis[1] + 1
        trim_EP2[n][2] = trim_EP2[n][2] - calibrate_axis[2] + 1
    return trim_EP2
"""

"""
def proj_xy_2d(path_of_sfam): # Old ver of proj
    target_sfam = open(path_of_sfam, 'rt', encoding = 'utf8')
    read_sfam = target_sfam.read()[:1500]
    target_sfam.close()
    read_sfam.find('define Lattice')
    for item in read_sfam.split("\n"):
        if 'define Lattice' in item:
            key_lattice = item.strip()
        else:
            pass
    x_latti = int(key_lattice.split(' ')[-3]) #436
    y_latti = int(key_lattice.split(' ')[-2]) #328
    z_latti = int(key_lattice.split(' ')[-1]) #74
    for item in read_sfam.split("\n"):
        if 'Data section follows' in item:
            idx_dsf = read_sfam.split("\n").index(item.strip())
        else:
            pass
    if read_sfam.split("\n")[idx_dsf+1][0] == '@':
        pre_matrix = read_sfam.split("\n")[idx_dsf+2:]
    else:
        pass
    pre_ar1 = np.array(list(map(int, list(filter(None, pre_matrix)))))
    peel_mx = pre_ar1.reshape(z_latti, y_latti, x_latti)
    main_mx = np.transpose(peel_mx, (2, 1, 0)) # (436, 328, 74)
    ijz_li = []
    proj_face_xy = np.zeros((x_latti, y_latti), dtype=int)
    for i in range(x_latti):
        for j in range(y_latti):
            for k in range(z_latti):
                ijz_li.append(main_mx[i][j][k])
            proj_face_xy[i][j] = max(ijz_li)
            ijz_li = []
    return proj_face_xy
"""

"""
def rdm_pick_slice_no(number, rdm_list): # Too slow
    record_name = rdm_list[number]
    count = number - rdm_list.index(rdm_list[number])
    return [record_name, count]
"""

"""
def blsoco(soma_2D_coordi_li): #blsoco stands for BLock SOma COordinate
	block_x_posi = soma_2D_coordi_li[0] // 100
	block_y_posi = soma_2D_coordi_li[1] // 100
	soma_block = (block_x_posi, block_y_posi)
	return soma_block
"""

"""
def blexhush(original_matrix): # blexhush stands for BLock EXtend HUndred SHape
    hundred_matrix = np.pad(original_matrix,
                            ((0,100 - np.shape(original_matrix)[0] % 100),(0,100 - np.shape(original_matrix)[1] % 100)),
                            'constant',
                            constant_values=(0))
    hundred_slice_li = [ [0 for _ in range(int(np.shape(hundred_matrix)[1]/100))] for _ in range(int(np.shape(hundred_matrix)[0]/100))]
    for hs_i in range(int(np.shape(hundred_matrix)[0] / 100)): # 500 / 100 = 5
        for hs_j in range(int(np.shape(hundred_matrix)[1] / 100)): # 400 / 100 = 4
            hundred_slice_li[hs_i][hs_j] = hundred_matrix[hs_i:hs_i+100,hs_j:hs_j+100]
    return list(hundred_slice_li)
"""

"""
 def projection_plot(matrix):
    proj_for_plot = matrix # to replace the original line above
    plt.figure(figsize=(10,10))
    plt.imshow(proj_for_plot, cmap = "gray",  aspect = 'equal')
    plt.show()
    return print('Shape of matrix projection:',np.shape(proj_for_plot))
"""
