# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:32:04 2017

@author: Benjamin
"""

import math
import numpy as np
from collections import Counter

# -----------------------------------------------------------------------------

"""
def common_elements(list1, list2):
    return list(set(list1) & set(list2))
"""




def find_duplicates(input_list):
    commons_found = False
    common_elems = [[]]*len(input_list)
    #print("input_list" + str(input_list))
    for idx1, elem1 in enumerate(input_list):
        for idx2, elem2 in enumerate(input_list):
            #print("elem1: " + str(elem1))
            #print("elem2: " + str(elem2))
            if not idx1 == idx2:
                #print("commons: " + str(list(set(elem1) & set(elem2))) )
                new_commons = list(set(elem1) & set(elem2))
                common_elems[idx1] = common_elems[idx1] + new_commons
                if not common_elems[idx1] == []:
                    commons_found = True
                    #print("commons_found: " + str(commons_found))
    
    unique_common_elems = []
    for bn in common_elems:
        np_bn = np.array(bn)
        np_bn = np.unique(np_bn)
        unique_common_elems.append(np_bn.tolist())
    # print(unique_common_elems)
    combined_bns = []
    for bn in unique_common_elems:
        combined_bns += bn
    # print(combined_bns)
    cnt = Counter(combined_bns)
    # print(cnt)
    duplicates = 0
    for entry in cnt:
        duplicates += cnt[entry] - 1
    
    return (duplicates, common_elems)


def remove_all(input_list, elem_to_remove):
    """
    Returns a copy of "input_list" in which all instances of "elem_to_remove"
    have been removed.
    """
    idx_to_keep = []
    for idx, elem in enumerate(input_list):
        if not elem == elem_to_remove:
            idx_to_keep.append(idx)
    # Build new list "return_list":
    return_list = []
    for idx in idx_to_keep:
        return_list = return_list + list([input_list[idx]])

    return return_list


def print_generations(bin_generations):
    for gen in bin_generations:
        line_str = ""
        for bn in gen:
            for sd in bn:
                line_str = line_str + str(sd)
            line_str = line_str + " -- "
        print(line_str)




