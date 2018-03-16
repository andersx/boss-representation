#!/usr/bin/env python

import numpy as np
import os
import qml

from qml.kernels import gaussian_kernel
from qml.kernels import laplacian_kernel
from qml.math import cho_solve

from qml.representations import get_slatm_mbtypes

from fboss import boss_kernel

from boss import generate_input

def get_energies(filename):
    """ Returns a dictionary with heats of formation for each xyz-file.
    """

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    energies = dict()

    for line in lines:
        tokens = line.split()

        xyz_name = tokens[0]
        hof = float(tokens[1])

        energies[xyz_name] = hof

    return energies

if __name__ == "__main__":

    test_dir = os.path.dirname(os.path.realpath(__file__))

    # Parse file containing PBE0/def2-TZVP heats of formation and xyz filenames
    data = get_energies(test_dir + "/hof_qm7.txt")

    # Generate a list of qml.Compound() objects
    mols = []


    keys = sorted(data.keys())
    
    # Shuffle molecules
    np.random.seed(666)
    np.random.shuffle(keys)

    n_test  = 500
    n_train = 1000
    
    n_total = n_test+n_train

    for xyz_file in keys[:n_total]:

        mol = qml.Compound(xyz=test_dir + "/qm7/" + xyz_file)
        mol.properties = data[xyz_file]
        mol.representation = generate_input(mol.nuclear_charges, mol.coordinates)
   
        mols.append(mol)


    # Make training and test sets

    training = mols[:n_train]
    test  = mols[-n_test:]

    # List of properties
    Y = np.array([mol.properties for mol in training])
    Ys = np.array([mol.properties for mol in test])
    
    X = np.array([mol.representation for mol in training])
    Xs = np.array([mol.representation for mol in test])

    print X.shape
    print Xs.shape


    # 1/r^n power law
    d_power = 4.0
   
    # Width for Gaussians in the spectrum
    d_width = 0.1

    # Gaussian-kernel width
    sigma = 50.0


    K  = boss_kernel(X, X, n_train, n_train, sigma, d_width, d_power)
    Ks = boss_kernel(X, Xs, n_train, n_test, sigma, d_width, d_power)
    
    print K

    K[np.diag_indices_from(K)] += 1e-8
    alpha = cho_solve(K,Y)

    # Calculate prediction kernel
    Yss = np.dot(Ks.transpose(), alpha)

    mae = np.mean(np.abs(Ys - Yss))

    print "boss"
    print mae
    
    for mol in training:
        mol.generate_bob()
    
    for mol in test:
        mol.generate_bob()

    X  = np.array([mol.representation for mol in training])
    Xs = np.array([mol.representation for mol in test])
   
    sigma = 100000.0
    llambda = 1e-10
    K  = laplacian_kernel(X, X, sigma)
    Ks = laplacian_kernel(X, Xs, sigma)
    
    K[np.diag_indices_from(K)] += llambda

    alpha = cho_solve(K,Y)

    # Calculate prediction kernel
    Yss = np.dot(Ks.transpose(), alpha)

    mae = np.mean(np.abs(Ys - Yss))

    print "bob"
    print mae
