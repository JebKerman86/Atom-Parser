# -*- coding: utf-8 -*-
"""
Einlesen von .XYZ Dateien
"""

import numpy as np
import json
import os
from pathlib import Path


CACHED_DATA_FOLDER_NAME = "cached_data"
INPUT_FOLDER_NAME = "input_files"
OUTPUT_FOLDER_NAME = "output_files"

# -----------------------------------------------------------------------------


def read_xyz_file(input_file_name):
    """
    Reads File "input_file_name".xyz, and returns two lists,
    specifying Atom types and Atom positions, respectively.
    """
    input_file_path = "./" + INPUT_FOLDER_NAME + "/" + input_file_name + ".xyz"
    file = open(input_file_path, 'r')
    max_file_lines = 1000

    num_atoms = int(file.readline())
    # file_comment = file.readline()
    atom_positions = []
    atom_types = []

    # skip comment section in file
    line = file.readline()

    iterations = 0
    line = file.readline()
    stripped_line = line.replace(" ", "").replace("\n", "")

    while stripped_line != "" and iterations < max_file_lines:

        line_string = line.split()
        atom_type_string = line_string[0]
        atom_types.append(atom_type_string)
        coord_string = line_string[1:]
        coord_xyz = [float(i) for i in coord_string]
        atom_positions.append(coord_xyz)

        iterations += 1
        line = file.readline()
        stripped_line = line.replace(" ", "").replace("\n", "")

    print("Length of Atom List: " + str(len(atom_positions)))

    if iterations >= max_file_lines:
        print("Input file too long!")

    if num_atoms != len(atom_positions):
        print("falsche Atomanzahl: num_atoms != Atoms.count")

    return (atom_types, atom_positions)


# -----------------------------------------------------------------------------


def read_transport_file():
    """
    Reads File "input_file_name".dat, and returns lists containing the atom
    indices of the device atoms, as well as the atom indices of
    the contact atoms. Also, a dictionary "interaction_distances" is generated,
    which spcifies the maximum interaction distance between each type of atom.
    """

    transport_file_path = "./" + INPUT_FOLDER_NAME + "/" + "transport.dat"
    file = open(transport_file_path, 'r')
    max_file_lines = 1000

    iterations = 0

    # IMPORTANT: In file, first atom has index is one, but in my program,
    # first atom has index is zero

    region_list = []  # List of regions, starting with device region

    line = file.readline()
    entries = line.split()
    

    iterations = 0
    while iterations < max_file_lines:

        region = list(range(int(entries[1]) - 1, int(entries[2])))
        region_list.append(region)

        line = file.readline()
        entries = line.split()
        iterations += 1
        
        if not("Device" in entries[0] or "Contact" in entries[0]):
            break

    interaction_distances = {}

    #line = file.readline()
    #stripped_line = line.replace(" ", "").replace("\n", "")
    #entries = line.split()

    # loop terminates at first empty line, or at end of file
    # (since readline() returns empty string at end of file)
    iterations = 0
    while iterations < max_file_lines:
        key = entries[0] + entries[1]
        interaction_distances[key] = float(entries[2])

        line = file.readline()
        entries = line.split()
        iterations += 1
        stripped_line = line.replace(" ", "").replace("\n", "")
        
        if stripped_line == '':
            break

    return (region_list, interaction_distances)

# -----------------------------------------------------------------------------


def write_bins(bin_generations, atom_positions, file_name):

    num_atoms = len(atom_positions)
    file_name = file_name + ".jmol"
    file_path = "./" + OUTPUT_FOLDER_NAME + "/" + file_name
    outfile = open(file_path, 'w')

    outfile.write(str(num_atoms) + "\n")
    outfile.write("This file shows the principal layers into which the molecule has been sorted." + "\n")
    
    chem_elements = ["H", "O"]
    
    #line_num = 3
    for gen_idx, generation in enumerate(bin_generations):
        for bn in generation:
            for atom_idx in bn:
                # print("atom_idx: " + str(atom_idx))
                line_str = chem_elements[gen_idx%2] + "    "
                for coord in atom_positions[int(atom_idx)]:
                    line_str = line_str + str(coord) + "    "
                outfile.write(line_str + "\n")
                #line_num += 1


# -----------------------------------------------------------------------------


def chache_data(input_file, data):
    """
    Caches matrices from function PrepInput() to lower execution time when
    program called repeatedly with the same input files.
    Converts Numpy arrays to lists.
    """
    jsonified_data = {}
    for key, array in data.items():
        converted_array = array
        if(str(type(array)) == "<class 'numpy.ndarray'>"):
            converted_array = array.tolist()
        jsonified_data[key] = converted_array

    input_file_name = input_file + ".json"
    file_path = "./" + CACHED_DATA_FOLDER_NAME + "/" + input_file_name
    with open(file_path, 'w') as outfile:
        json.dump(jsonified_data, outfile)


# -----------------------------------------------------------------------------


def load_data(input_file_name):
    """
    Loads cached data for reuse. Converts Lists to numpy arrays.
    """

    jsonified_data = {}
    data = {}
    file_name = input_file_name + ".json"
    file_path = "./" + CACHED_DATA_FOLDER_NAME + "/" + file_name
    my_file = Path(file_path)

    if my_file.is_file():

        with open(str(file_path), 'r') as file:
            jsonified_data = json.load(file)

            for key, array in jsonified_data.items():
                converted_array = np.array(array)
                data[key] = converted_array

            return data
    
    print("No File Found")
    return "No File Found"
