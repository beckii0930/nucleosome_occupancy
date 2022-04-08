import numpy as np;
import time
import scipy.io
import sys
import math

def readInputAsArray(fileName):
    with open(fileName, 'r') as myfile:
        data = myfile.readlines()

    # Strip newline
    for i in range(0, len(data)):
        data[i] = data[i].rstrip()
    # print(data)
    return data

def my_print(name, val):
    print(f'{name}: {val}')
    return

# test=np.array([10])
# my_print('test', test)

def list_str_to_float(arr):
    out = []
    for i in arr:
        out+=[float(i)]
    return out

def append_list_to_np_arr(lst, np_arr):
    bs_lst=[]
    for i in range(len(lst)):
#         bs_lst+= [[lst[i], lst[i]]]
        bs_lst+= [[lst[i]]]
    
    np_lst=np.array([bs_lst])
#     my_print('np_lst.shape', np_lst.shape)

    if np_arr.shape[0] == 0:
        np_arr = np_lst 
    else:
        np_arr = np.concatenate((np_arr, np_lst), axis=0) 
#     my_print('np_lst.shape', np_arr.shape)
    return np_arr

def append_bplist_to_np_arr(lst, np_arr):
    bp_lst = []
    for i in range(1, len(lst)):
        bp_lst+= [[lst[i-1], lst[i]]]
    np_bp=np.array([bp_lst])
    
#     my_print('np_bp', np_bp.shape)
#     my_print('np_arr', np_arr.shape)
    if np_arr.shape[0] == 0:
        np_arr = np_bp
    else:
        np_arr = np.concatenate((np_arr, np_bp), axis=0) 
#     print(np_arr.shape)
    return np_arr


def makeSingleShapeArr(seq_file, total_sections, section):
#     print(">>>>>>>>> In making makeSingleShapeArr ")
    np_all_shape_arr = np.array([])
    
    Seqs = readInputAsArray(seq_file)
    seq_count=1
    bp_shape=False
    total_lines = 0;
    total_sections = int(total_sections);
    section = int(section);
    
    for line in Seqs:
        l = line.split()
        if l[1] != '3': # Only consider flanking region of size 3
            continue;
        total_lines += 1;
        
    print(f"Total # of regions is: {total_sections}\n");
    print(f"Current regions is: {section}\n");
    print(f"Total # of lines is: {total_lines}\n");
    section_length = math.floor(total_lines / total_sections);
    start_line = (section-1) * section_length;
    end_line = section * section_length-1;

    if (end_line > total_lines):
        end_line = total_lines-1;

    line_count = 0;
    for line in Seqs:
        l = line.split()
        if l[1] != '3': # Only consider flanking region of size 3
            continue;
        if(line_count < start_line):
            line_count += 1;
            continue;
        if(line_count > end_line):
            break;
        line_count += 1;
        
        curr_seq_shape_list = list_str_to_float(l[2:]) # get the shape values only
        if len(l[0]) > len(curr_seq_shape_list):
#             print("base step parameter")
            bp_shape = False
        else:
#             print("base pair parameter")
            bp_shape = True
        
        ## 1. for base pair shapes, there are 147 vals for each seq 
        if bp_shape: 
            if len(curr_seq_shape_list) < 147:
                while len(curr_seq_shape_list) < 147: # if the sequence is shorter than 147
                    curr_seq_shape_list+=[0]
                    
                np_all_shape_arr = append_bplist_to_np_arr(curr_seq_shape_list, np_all_shape_arr)
            else:
                start = 0
                 # is a base pair shape, each 147bp seq has 147 shape vals, we use 147-2+1
                while (start < len(curr_seq_shape_list)-146):
                    curr_start = start;
                    curr_end = 146 + start;
                    sub_seq_shape_list = curr_seq_shape_list[curr_start: curr_end+1];
                    start+=1
                    np_all_shape_arr = append_bplist_to_np_arr(sub_seq_shape_list, np_all_shape_arr)

        ## 2. for base step shapes, there are 146 vals for each seq 
        else:
            if len(curr_seq_shape_list) < 146:
                while len(curr_seq_shape_list) < 146: # if the sequence is shorter than 146
                    curr_seq_shape_list+=[0]
                    
                np_all_shape_arr = append_list_to_np_arr(curr_seq_shape_list, np_all_shape_arr)
            else:
#                 print('>>>>>> Found longer sequences')
                start = 0
                # is a base step shape, each 147bp seq has 146 shape vals, we use 146-2+1
                while (start < len(curr_seq_shape_list)-145):
                    curr_start = start;
                    curr_end = 145 + start;
                    sub_seq_shape_list = curr_seq_shape_list[curr_start: curr_end+1];
                    start+=1
                    np_all_shape_arr = append_list_to_np_arr(sub_seq_shape_list, np_all_shape_arr)
#         my_print('np_all_shape_arr', np_all_shape_arr.shape)
        seq_count+=1
#     my_print('np_all_shape_arr.shape', np_all_shape_arr.shape)
    return np_all_shape_arr;


###### Now merge the single shape arrs to a whole np array
def concatAllShapeArr(seq_file, All_Shapes, total_sections, section):
    np_shape_E = np.array([]) # enriched seqs output

    for i in range(len(All_Shapes)):
#         my_print('\ncurr shape is ', All_Shapes[i])
        np_all_shape_arr = makeSingleShapeArr(seq_file + All_Shapes[i] + '.txt', total_sections, section)
        if np_shape_E.shape[0] == 0:
            np_shape_E = np_all_shape_arr
        else:
            ## merge exp1: consider the base step as separate features
            np_shape_E = np.concatenate([np_shape_E, np_all_shape_arr], axis=2) #(218, 146, 3)
            
#             np_shape_E = np.stack([np_shape_E, np_all_shape_arr], axis=2) #(218, 146, 2, 2)
    return np_shape_E

################################################## Main ##################################################
# seq_file='/home/yibei/Downloads/yibei_predictions/'
seq_file='/project/rohs_108/yibeijia/data/yibei_predictions/'
All_Shapes=['Buckle-FL', 'Buckle', 'EP', 'HelT-FL', 'HelT', 'MGW-FL', 'MGW', 
                'Opening-FL', 'Opening', 'ProT-FL', 'ProT', 'Rise-FL', 'Rise', 'Roll-FL',
                'Roll', 'Shear-FL', 'Shear', 'Shift-FL', 'Shift', 'Slide-FL', 'Slide', 
                'Stagger-FL', 'Stagger', 'Stretch-FL', 'Stretch', 'Tilt-FL', 'Tilt']

# tb_file='/home/yibei/Projects/nucleosome_occupancy/pentamer_table.txt'
# seq_file='/Users/yibeijia/Downloads/yibei_predictions/'
np_shape_E = np.array([])
np_shape_D = np.array([])

total_sections=sys.argv[1]
section = sys.argv[2]
# All_Shapes=['Stretch']
print("Getting Enriched Seqs")


np_shape_E = concatAllShapeArr(seq_file+'enriched_', All_Shapes, total_sections, section) # depleted seqs output

# All_Shapes=['Roll', 'Stretch']
print("Getting Depleted Seqs")
np_shape_D = concatAllShapeArr(seq_file+'depleted_', All_Shapes, total_sections, section) # enriched seqs output

my_print('np_shape_E', np_shape_E.shape) # (NumSeq, 147, 39)
my_print('np_shape_D', np_shape_D.shape) # (NumSeq, 147, 39)

data = {"EnrichedData": np_shape_E, "DepletedData": np_shape_D};
# mat_filename = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/sampleSeqs.mat';
# mat_filename = '/home/yibei/Projects/data/train_test_data/ShapeDNNSeqs_' + str(section) + '.mat';
mat_filename = '/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/ShapeDNNSeqs_' + str(section) + '.mat';

scipy.io.savemat(mat_filename, data,  do_compression=True);
