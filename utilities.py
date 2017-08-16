# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:32:04 2017

@author: Benjamin
"""

import math


# -----------------------------------------------------------------------------

"""
def common_elements(list1, list2):
    return list(set(list1) & set(list2))
"""

def common_elements(input_list):
    commons_found = False
    common_elems = [[]]*len(input_list)
    #print("input_list" + str(input_list))
    for idx1, elem1 in enumerate(input_list):
        for idx2, elem2 in enumerate(input_list):
            #print("elem1: " + str(elem1))
            #print("elem2: " + str(elem2))
            if not idx1 == idx2:
                print("commons: " + str(list(set(elem1) & set(elem2))) )
                common_elems[idx1] = common_elems[idx1] + list(set(elem1) & set(elem2))
                if not common_elems[idx1] == []:
                    commons_found = True
                    print("commons_found: " + str(commons_found))

    return (commons_found, common_elems)


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
