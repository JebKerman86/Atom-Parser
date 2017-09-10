# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 14:32:04 2017

@author: Benjamin
"""

from bin_sort import bins_are_neighbours

# -----------------------------------------------------------------------------


def print_generations(bin_generations):
    for gen_idx, gen in enumerate(bin_generations):
        # print("gen_idx: " + str(gen_idx))
        # print(gen)
        line_str = ""
        for bn in gen:
            # print("bn: " + str(bn))
            line_str = line_str + str([x+1 for x in bn])
            line_str = line_str + " -- "
        print(line_str)


def print_final_chain(final_chain):
    line_str = ""
    for bn in final_chain:
        line_str = line_str + str([x+1 for x in bn])
        line_str = line_str + "\n"

    print(line_str)


def count_atoms(chains):
    num_atoms = 0
    sorted_atom_idxs = []
    for gen in chains:
        for bn in gen:
            for atom_idx in bn:
                if atom_idx not in sorted_atom_idxs:
                    sorted_atom_idxs.append(atom_idx)
                    num_atoms += 1
                else:
                    print("In count_atoms: duplicate atom found.")
    return num_atoms


def test_solution(final_chain, interact_mtrx):

    # Check whether all bins have two neighbour bins
    num_neighbours = [0]*len(final_chain)
    # print("num_neighbours: " + str(num_neighbours))
    for gen_idx1, gen1 in enumerate(final_chain):
        for gen_idx2, gen2 in enumerate(final_chain):
            if not gen_idx2 == gen_idx1:
                if bins_are_neighbours(gen1, gen2, interact_mtrx):
                    num_neighbours[gen_idx1] += 1
    print("num_neighbours: " + str(num_neighbours))
    for num_n in num_neighbours:
        if num_n > 2:
            print("At least one bin has more than two neighbours.")
            print("---> BAD SOLUTION <---")
            return False

    # Check whether the atom_idx from all_atom_idxs exist exactly once
    # in the final_chain, and no other atom_idx are present in the final_chain.
    all_atom_idxs = list(range(len(interact_mtrx[:, 0])))

    for gen in final_chain:
        for atom_idx in gen:
            if atom_idx in all_atom_idxs:
                all_atom_idxs.remove(atom_idx)
            else:
                print("Either duplicate atom_idx, \
                      or atom_idx not part of original molecule found:")
                print("atom_idx: " + str(atom_idx+1))
                print("---> BAD SOLUTION <---")
                return False
    if not len(all_atom_idxs) == 0:
        print("At least one atom_idx from original molecule was \
              not sorted into final_chain:")
        print("unsorted atom_idx: " + str([x+1 for x in all_atom_idxs]))
        print("---> BAD SOLUTION <---")
        return False

    print("Solution check passed!")
    return True
