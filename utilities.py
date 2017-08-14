# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:32:04 2017

@author: Benjamin
"""

import math


# -----------------------------------------------------------------------------


def print_list(input_list):
    """
    Prints a list! (duh)
    """
    print("--------------------------------------------------------")
    for elem in input_list:
        print(elem)
    print("--------------------------------------------------------")


# -----------------------------------------------------------------------------


def print_matrix(input_matrix):
    """
    Prints a matrix! (duh)
    """
    print("------------------------------------------------------------------")
    for line in input_matrix:
        print_line = ""
        for number in line:
            print_line = print_line + "{:5.2f}".format(number) + "  |  "
        print(print_line)
    print("------------------------------------------------------------------")


# -----------------------------------------------------------------------------


def vec3_dist(vec1, vec2):
    """
    Returns distance between Vectors (list of numbers) vec1 and vec2
    """
    sum_sqr_dist = 0

    for dimension in range(0, len(vec1)):
        sum_sqr_dist += math.pow(vec1[dimension] - vec2[dimension], 2)

    return math.sqrt(sum_sqr_dist)


# -----------------------------------------------------------------------------


# Numpy: sort, argsort

# Returns a list ("found_List") of the indices of the smallest element,
# ie if only one smallest element => len(found_list) == 1).
# Indices in found_list in same order as they oeccured in th input_list
def smallest_element_indices(input_list, indices_of_candidate_elements):
    """
    In input_list, find indices of smallest elements (can be multiple if
    there are duplicate elements). Return list of these indices.
    indices_of_candidate_elements is a list of indices that are open for
    consideration, all atom indices NOT on that list are ignored.
    """
    found_list = []

    if not input_list:
        return found_list

    elem, idx = min((input_list[i], i) for i in indices_of_candidate_elements)
    elem_old = elem

    # while loop runs as long as there is a duplicate of the smallest element
    while elem_old == elem:
        found_list.append(idx)
        indices_of_candidate_elements.remove(idx)
        elem_old = elem
        if not indices_of_candidate_elements:  # Exit loop if list is empty
            break
        elem, idx =  \
            min((input_list[i], i) for i in indices_of_candidate_elements)

    return found_list


# -----------------------------------------------------------------------------


def order_index_list(input_list):
    """
    Take a list of distances in "input_list", and return indices corresponding
    to those distances in the list "result", ordered by ascending distances
    """

    ordered_list = []
    indices_of_candidate_elements = list(range(len(input_list)))

    while indices_of_candidate_elements:
        ordered_list.append(
            smallest_element_indices(input_list,
                                     indices_of_candidate_elements))

    result = []
    # unpack the list
    for _list in ordered_list:
        for elem in _list:
            result.append(elem)

    return result





def np_order(input):
    return input








