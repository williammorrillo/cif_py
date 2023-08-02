"""
Utility functions for the package.
"""

import numpy as np
import re

def lattice_matrix(data):
    """
    Returns the lattice matrix from the cif data.
    """
    a = float(data['_cell_length_a'][0])
    b = float(data['_cell_length_b'][0])
    c = float(data['_cell_length_c'][0])
    alpha = float(data['_cell_angle_alpha'][0])
    beta = float(data['_cell_angle_beta'][0])
    gamma = float(data['_cell_angle_gamma'][0])

    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    gamma = np.deg2rad(gamma)

    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)

    a1 = [a, 0, 0]
    a2 = [b*cos_gamma, b*sin_gamma, 0]
    a3 = [c*cos_beta, c*(cos_alpha-cos_beta*cos_gamma)/sin_gamma, c*np.sqrt(1+2*cos_alpha*cos_beta*cos_gamma-cos_alpha**2-cos_beta**2-cos_gamma**2)/sin_gamma]

    return np.array([a1, a2, a3])


def remove_occupancy(list):
    """
    Removes the () from a number in a list.
    """
    return [x.split('(')[0] for x in list]


def get_labels(data):
    """
    returns the atoms and labels from the cif data.
    """
    atom_list = data['_atom_site_label']
    atom_list = [re.split(r'(\d+)', x)[0] for x in atom_list]
    return atom_list


def sort_atoms(coords):
    coords = coords[coords[:, 0].argsort()]
    atoms = []
    natoms = []
    atom_list = coords[:, 0]
    for atom in atom_list:
        if atom not in atoms:
            atoms.append(atom)
            natoms.append(1)
        else:
            natoms[atoms.index(atom)] += 1

    atoms = [(atom, natom) for atom, natom in zip(atoms, natoms)]
    atoms = sorted(atoms, key=lambda x: x[0])

    natoms = [str(x[1]) for x in atoms]
    atoms = [x[0] for x in atoms]
    return atoms, natoms, coords


def get_sym_ops(data):
    """
    Performs symmetry operations on the coords
    """
    sym_ops = None

    try:
        sym_ops = data['_symmetry_equiv_pos_as_xyz']
    except:
        pass

    try:
        sym_ops = data['_space_group_symop_operation_xyz']
    except:
        pass

    if sym_ops is not None:

        replace_str_int(sym_ops)

        return sym_ops


def replace_str_int(sym_ops):
    for i in range(len(sym_ops)):
        for j in range(3):
            sym_ops[i, j] = sym_ops[i, j].replace("'", "")
            sym_ops[i, j] = sym_ops[i, j].replace(",", "")
            sym_ops[i, j] = sym_ops[i, j].replace("x", "1")
            sym_ops[i, j] = sym_ops[i, j].replace("y", "1")
            sym_ops[i, j] = sym_ops[i, j].replace("z", "1")
            sym_ops[i, j] = sym_ops[i, j].replace("1/2", "0.5")
            sym_ops[i, j] = sym_ops[i, j].replace("1/3", "0.3333333333333333")
            sym_ops[i, j] = sym_ops[i, j].replace("2/3", "0.6666666666666666")
            sym_ops[i, j] = sym_ops[i, j].replace("1/4", "0.25")
            sym_ops[i, j] = sym_ops[i, j].replace("3/4", "0.75")
            sym_ops[i, j] = sym_ops[i, j].replace("1/5", "0.2")
            sym_ops[i, j] = sym_ops[i, j].replace("2/5", "0.4")
            sym_ops[i, j] = sym_ops[i, j].replace("3/5", "0.6")
            sym_ops[i, j] = sym_ops[i, j].replace("4/5", "0.8")
            sym_ops[i, j] = sym_ops[i, j].replace("1/6", "0.16666666666666666")
            sym_ops[i, j] = sym_ops[i, j].replace("5/6", "0.8333333333333333")
            sym_ops[i, j] = sym_ops[i, j].replace("1/7", "0.14285714285714285")
            sym_ops[i, j] = sym_ops[i, j].replace("2/7", "0.2857142857142857")
            sym_ops[i, j] = sym_ops[i, j].replace("3/7", "0.42857142857142855")
            sym_ops[i, j] = sym_ops[i, j].replace("4/7", "0.5714285714285714")
            sym_ops[i, j] = sym_ops[i, j].replace("5/7", "0.7142857142857143")
            sym_ops[i, j] = sym_ops[i, j].replace("6/7", "0.8571428571428571")
            sym_ops[i, j] = sym_ops[i, j].replace("1/8", "0.125")
            sym_ops[i, j] = sym_ops[i, j].replace("3/8", "0.375")
            sym_ops[i, j] = sym_ops[i, j].replace("5/8", "0.625")
            sym_ops[i, j] = sym_ops[i, j].replace("7/8", "0.875")

            if sym_ops[i, j].startswith("-"):
                # split by position 2
                sym_ops[i, j] = (sym_ops[i, j][0:2], sym_ops[i, j][3:])
            else:
                # split by position 1
                sym_ops[i, j] = (sym_ops[i, j][0:1], sym_ops[i, j][2:])

            if sym_ops[i, j][1] == '':
                sym_ops[i, j] = (sym_ops[i, j][0], 0)

    return sym_ops


def apply_sym_ops(sym_ops, coords):

    new_coords = np.array([[x[0],
                            float(x[1]) * float(y[0][0]) + float(y[0][1]),
                            float(x[2]) * float(y[1][0]) + float(y[1][1]),
                            float(x[3]) * float(y[2][0]) + float(y[2][1])]
                           for x in coords
                           for y in sym_ops])

    new_coords = np.unique(new_coords, axis=0)

    return new_coords
