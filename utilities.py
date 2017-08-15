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
    any_commons = False
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
                    any_commons = True
                    print("any_commons: " + str(any_commons))

    return (any_commons, common_elems)

