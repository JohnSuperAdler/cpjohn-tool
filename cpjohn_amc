##### cpjohn_amc.py
#####
##### "amc" stands for AM file configureration.
#####
##### It's a manual module made up for Analysis of Drosophila brain neuron
##### 3-D matrix data. This py module file contents functions below: 


import numpy as np
import matplotlib.pyplot as plt
import os
import re


##### Most used

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

def mx_preview(array):
    #plt.figure(figsize=(8, 8))
    plt.imshow(array, aspect = 'equal', interpolation='nearest')
    plt.title('Preview of(x, y) = {}'.format(np.shape(array)),
              fontsize = 20, y=1.05)
    plt.xlabel('Y', fontsize = 20)
    plt.ylabel('X', fontsize = 20)
    plt.colorbar(fraction=0.04, pad=0.1)
    return plt.show()

def happy_time(start,stop):
    process_time = round(stop - start)
    ss = process_time % 60
    mm = process_time // 60 % 60
    hh = process_time // 3600
    duration = "Process time == {}s == {}H {}m {}s".format(process_time, hh, mm, ss)
    return duration

##### Actually not used very often in these months

def soma_index(path_of_branches_txt):
    # Part of loading soma position (proto) in branch.txt
    brtxt_filename = str(path_of_branches_txt)
    target_bptxt = open(brtxt_filename, 'rt', encoding = 'utf8')
    read_bptxt = target_bptxt.read()[:1000]
    target_bptxt.close()
    trim_soma = re.findall('--*.\n[0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+ [0-9]+\n--*.\n',read_bptxt)
    soma_index_proto = list(map(int, trim_soma[0].split('\n')[1].split(' ')[-5::-1]))
    return soma_index_proto

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

def rdm_pick_slice_no(number, rdm_list):
    record_name = rdm_list[number]
    count = 0
    while rdm_list[number] == rdm_list[number-1]:
        count += 1
        number += -1
    return [record_name, count]

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
