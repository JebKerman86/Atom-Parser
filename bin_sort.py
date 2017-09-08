# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 15:11:13 2017

@author: Benjamin
"""

import numpy as np
from copy import deepcopy

from utilities import remove_all, print_generations

def get_contact_bins(device, contacts, interact_mtrx):

    contact_starter_atom_list = []

    for contact in contacts:
        contact_edge_list = []
        for index_contact in contact:
            for index_device in device:
                if interact_mtrx[index_contact, index_device]:
                    if not (index_contact in contact_edge_list):
                        contact_edge_list.append(index_contact)
        contact_starter_atom_list.append(np.array(contact_edge_list))

    return contact_starter_atom_list



def get_next_bins(curr_bins, prev_bins, interact_mtrx):
    """
    List of bins. Each element contains the bins corresponding to a contact.
    """
    next_bins = []
    num_atoms = np.size(interact_mtrx, axis=0)
    atoms = np.array(list(range(num_atoms)))

    for curr_bin in curr_bins:
        bin_candidates = np.array([])
        for atom_idx in curr_bin:
            # print(atom_idx)
            # print(interact_mtrx[atom_idx, :])
            bin_add_candidates = atoms[interact_mtrx[atom_idx, :]]
            for prev_bin in prev_bins:
                #print("prev_bin" + str(prev_bin))
                bin_add_candidates = [x for x in bin_add_candidates if x not in prev_bin]
                bin_add_candidates = np.array(bin_add_candidates)
            bin_candidates = np.append(bin_candidates, bin_add_candidates)
            bin_candidates = [int(x) for x in bin_candidates]
            # print("bin_add_candidates" + str(bin_add_candidates))
            # print("bin_candidates" + str(bin_candidates))
        bin_atoms = np.unique(bin_candidates)
        #print("bin_atoms: " + str(bin_atoms))
        next_bins.append(bin_atoms)
    
    return(next_bins)


def remove_common_elems(prev_gen, curr_gen, common_elems_by_chain):
    """
    Takes the bins at tips of each contact chain and ensures that there are no
    duplicate atoms in these bins. The atom in the contact chain with
    the smallest index is kept.
    """
    atoms_to_delete = []
    for chain_idx, common_atoms in enumerate(common_elems_by_chain):
        for idx, atom_idx in enumerate(common_atoms):
            if not atom_idx in atoms_to_delete:
                atoms_to_delete.append(atom_idx)
                print("add atom_idx: " + str(atom_idx))
                prev_gen[chain_idx] = \
                     np.append(prev_gen[chain_idx], np.array([atom_idx]) )

        for atom_idx in atoms_to_delete:
            # print("elems_to_delete: " + str(elems_to_delete))
            curr_gen[chain_idx] =  \
                remove_all(curr_gen[chain_idx], atom_idx)



def bins_are_neighbours(bin1, bin2, interact_mtrx):
    """
    Check if bin1 and bin2 contain interacting atoms, or have atoms in common.
    """

    for atom_idx1 in bin1:
        for atom_idx2 in bin2:
            #if atom_idx2 > atom_idx1:   (????????)
            if interact_mtrx[atom_idx1, atom_idx2]:
                return True
    return False




def check_for_collisions(bin_generations, gen_idx1, gen_idx2, interact_mtrx):
    # print("Entering check_for_collisions...")
    # print("gen_idx1 = " + str(gen_idx1))
    # print("gen_idx2 = " + str(gen_idx2))
    
    gen1 = bin_generations[gen_idx1]
    gen2 = bin_generations[gen_idx2]
    # print("gen1: " + str(gen1))
    # print("gen2" + str(gen2))
    
    collisions = []
    for bin_idx1, bin1 in enumerate(gen1):
        for bin_idx2, bin2 in enumerate(gen2):
            if bin_idx2 > bin_idx1:
                col = bins_are_neighbours(bin1, bin2,interact_mtrx)
                if col:
                    print("collision detected: ")
                    print("bin1" + str(bin1))
                    print("bin2" + str(bin2))
                    collisions.append([(gen_idx1, gen_idx2),(bin_idx1, bin_idx2)])
                    
    # print("Exiting check_for_collisions.")
    if len(collisions) == 0:
        return (False, [(-1,-1)])
    else:
        return (True, collisions)



def check_generation_for_collisions(bin_generations, gen_idx, interact_mtrx):
    
    curr_generation = bin_generations[gen_idx]
    if gen_idx-1 > 0:
        prev_generation = bin_generations[gen_idx-1]

    # Return list of collision tuples, since there can be multiple collisions
    # in a generation
    collisions = []
    for bin_idx1, bin1 in enumerate(curr_generation):
        for bin_idx2, bin2 in enumerate(curr_generation):
            if bin_idx2 > bin_idx1:
                if (not bin1 == []) and (not bin2 == []):
                    col = bins_are_neighbours(bin1, bin2,interact_mtrx)
                if bin1 == []:
                    print("check prev gen")
                    col = bins_are_neighbours(prev_generation[bin_idx1], bin2,interact_mtrx)
                if bin2 == []:
                    print("check prev gen")
                    col = bins_are_neighbours(bin1, prev_generation[bin_idx2],interact_mtrx)
                if col:
                    print("collision detected: ")
                    print("bin1" + str(bin1))
                    print("bin2" + str(bin2))
                    collisions.append((bin_idx1, bin_idx2))
                    
    if len(collisions) == 0:
        return (False, [(-1,-1)])
    else:
        return (True, collisions)


def find_all_collisions(bin_generations, interact_mtrx):
    
    print("Entering find_all_collisions... ")
    # col_list is list of tuples:
    # [ (gen_idx1, [(chain_idx1,chain_idx2), (chain_idx1,chain_idx2)]),
    #   (gen_idx2, [(chain_idx1, chain_idx2)]) ]
    
    intergen_col_list = []
    for gen_idx1, gen1 in enumerate(bin_generations):
        for gen_idx2, gen2 in enumerate(bin_generations):
            if gen_idx2 >= gen_idx1:
                col = check_for_collisions(bin_generations, gen_idx1, gen_idx2, interact_mtrx)
                if col[0]:
                    intergen_col_list.append(col)
    
    
    print("intergen_col_list: ")
    print(intergen_col_list)
    
    col_list = []
    num_chains = len(bin_generations[0])
    all_chain_idxs = list(range(num_chains))
    # "-1" in uncollided_chain_idxs is a dummy value that allows while loop to begin
    uncollided_chain_idxs = [-1]
    iterations = 1
    
    
    """
    DO I NEED THIS WHILE LOOP ?????!!?
    """
    """
    while not uncollided_chain_idxs == [] and iterations < 10:
        iterations += 1
        col_list = []
    """

    for gen_idx in range(len(bin_generations)):
        col = check_generation_for_collisions(bin_generations, gen_idx, interact_mtrx)
        # print("gen_idx = " + str(gen_idx) + "  /  col = " + str(col))
        if col[0]:
            col_list.append((gen_idx, col[1]))
    """
        # Figure out whether any chains have NOT collided with eachother
        collided_chain_idxs = []
        print("col_list" + str(col_list))
        for col in col_list:
            col_tuple = col[1][0]
            print("col_tuple: " + str(col_tuple))
            collided_chain_idxs.append(col_tuple[0])
            collided_chain_idxs.append(col_tuple[1])
        print("collided_chain_idxs: " + str(collided_chain_idxs))
    
        uncollided_chain_idxs = []
        for chain_idx in all_chain_idxs:
            if not chain_idx in collided_chain_idxs:
                uncollided_chain_idxs.append(chain_idx)
        print("uncollided_chain_idxs: " + str(uncollided_chain_idxs))
        
        # For all uncollided chains, merge the last bin into the second to
        # last bin.
        for chain_idx in uncollided_chain_idxs:
            merge_chain_tip(bin_generations, chain_idx)

        print("In find_all_collisions: bin_generations: ")
        print_generations(bin_generations)
    """
    print("Exiting find_all_collisions... ")
    return col_list


def find_chain_tip_idx(bin_generations, chain_idx):
    for gen_idx, gen in enumerate(bin_generations):
        bn = gen[chain_idx]
        if len(bn) == 0:
            last_full_bin_gen_idx = gen_idx-1
            return last_full_bin_gen_idx
        

def merge_chain_tip(bin_generations, chain_idx):
    chain_tip_gen_idx = find_chain_tip_idx(bin_generations, chain_idx)

    for atom_idx in bin_generations[chain_tip_gen_idx][chain_idx]:
        bin_generations[chain_tip_gen_idx-1][chain_idx] = \
            np.append(bin_generations[chain_tip_gen_idx-1][chain_idx], np.array([atom_idx]))
    bin_generations[chain_tip_gen_idx][chain_idx] = np.array([])
    
    new_tip_gen_idx = chain_tip_gen_idx - 1
    return new_tip_gen_idx


def check_collision(collision_list, chain_idx, chain_tip_gen_idx):

    for collision in collision_list:
        if chain_tip_gen_idx == collision[0]:
            for col_tuple in collision[1]:
                if chain_idx in col_tuple:
                    return True

    return False


def trim_until_collision(bin_generations, collision_list, chain_idx):

    has_collided = False
    tip_gen_idx = find_chain_tip_idx(bin_generations, chain_idx)
    while tip_gen_idx > 0:
        has_collided = check_collision(collision_list, chain_idx, tip_gen_idx)
        if has_collided == True:
            break
        merge_chain_tip(bin_generations, chain_idx)
        tip_gen_idx = find_chain_tip_idx(bin_generations, chain_idx)





def merge(chains, col_gen_idx, chain1_idx, chain2_idx):
    
    merged_chains = []
    for gen in chains:
        merged_chains.append(deepcopy(gen))

    # merged_chains.append(col_gen)

    for gen_idx in range(len(merged_chains)):
        # print("gen_idx " + str(gen_idx))
        merged_bin = np.concatenate(
                (deepcopy(merged_chains[gen_idx][chain1_idx]),
                 deepcopy(merged_chains[gen_idx][chain2_idx])),
                 axis=0)
        merged_bin = np.unique(np.array(merged_bin))
        # merged_bin = [int(x) for x in merged_bin]
        print("merged_bin: ")
        print([x+1 for x in merged_bin])
        merged_chains[gen_idx][chain1_idx] = merged_bin
        merged_chains[gen_idx][chain2_idx] = np.array([])

    return merged_chains


def merge_chain(merged_bin_generations, bin_generations, collision_list, col_gen_idx, col_tuple):

    # Handle cases when more than two chains collide at same time at same place
    # col_tuple[0] always contains the smaller chain index
    # print("col_tuple")
    # print(col_tuple)
    # print("merged_bin_generations[0][2]: ")
    # print(merged_bin_generations[0][2] + merged_bin_generations[0][0])
    
    chain1_idx = col_tuple[0]
    chain2_idx = col_tuple[1]

    """
    Check whether any of the two collided chains have uncollided tips.
    If so, trim the chains down until the tip is involved in a collision
    """
    trim_until_collision(merged_bin_generations, collision_list, chain1_idx)
    trim_until_collision(merged_bin_generations, collision_list, chain2_idx)
    
    print("In merge_chain: bin_generations: ")
    print_generations(merged_bin_generations)

    
    for gen_idx in range(col_gen_idx+1):
        # print("gen_idx " + str(gen_idx))
        merged_bin = np.concatenate(
                (deepcopy(merged_bin_generations[gen_idx][chain1_idx]),
                 deepcopy(merged_bin_generations[gen_idx][chain2_idx])),
                 axis=0)
        merged_bin = [int(x) for x in merged_bin]
        print("merged_bin: ")
        print(merged_bin)
        merged_bin_generations[gen_idx][chain1_idx] = merged_bin
        merged_bin_generations[gen_idx][chain2_idx] = np.array([])
        
    # print("merged_bin_generations: ")
    # print(merged_bin_generations)



def glue_chains(chain1, chain2):
    """
    Glue dangling ends of two chains together:
    [0 --> chain1  -->  -1]  +  [-1  -->  chain2  -->  0]
    """
    
    new_chain = []
    last_bin = []

    for gen_idx, bn in enumerate(chain1):
        new_chain.append(deepcopy(bn))
        if gen_idx+1 == len(chain1):
            last_bin = new_chain[-1]
    print(last_bin)
    for ii in range(1, len(chain2)+1):
        idx = len(chain2) - ii
        bn = deepcopy(chain2[idx])
        bn = [x for x in bn if x not in last_bin]
        if len(bn) == 0:
            continue
        new_chain.append(bn)
        
    return new_chain



"""
def find_contacts_to_keep(bin_generations):
    
    Find the two contacts that should be left over after merging.
    These are the two contacts that have the most bins (principal layers)
    between them.
    
    num_ctcts = len(bin_generations[0])
    for cntct_idx1 in num_ctcts:
        for cntc_idx2 in num_ctcts:
            if cntc_idx2 > cntct_idx1:
                
"""


