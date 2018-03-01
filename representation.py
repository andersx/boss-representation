#!/usr/bin/env python

import numpy as np
import os
import qml

from qml.kernels import gaussian_kernel
from qml.math import cho_solve

from qml.representations import get_slatm_mbtypes

MAX_ID = 16
MAX_ONE_BODY = 16
MAX_TWO_BODY = 200
MAX_THREE_BODY = 200


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


def generate_one_body(coordinates, nuclear_charge):


    one_body = np.zeros((MAX_ID, MAX_ONE_BODY))
    one_body_counts = np.zeros(MAX_ID, dtype = np.int32)

    for i, nuc in enumerate(nuclear_charge):

        q = nuc - 1

        # print i, q

        one_body[q,one_body_counts[q]] += 1
        one_body_counts[q] += 1

    return one_body, one_body_counts

def generate_two_body(coordinates, nuclear_charge):
    
    n = 0

    indexes = dict()

    for i in range(1, MAX_ID +1):
        for j in range(i, MAX_ID +1):

            # print n, i,j

            key = (i,j)

            indexes[key] = n

            n += 1

    two_body_distance = np.zeros((n, MAX_TWO_BODY))
    two_body_scaling = np.zeros((n, MAX_TWO_BODY))
    two_body_counts = np.zeros(n, dtype = np.int32)


    for i, nuc1 in enumerate(nuclear_charge):

        q1 = nuc1 - 1
        for j, nuc2 in enumerate(nuclear_charge):
   
            if (i <= j): continue

            q2 = nuc2 - 1

            key = (min(nuc1, nuc2), max(nuc1, nuc2))
            n = indexes[key]

            dist = distance(coordinates[i], coordinates[j])
            scaling = dist**-4
            
            # print i, j, key, scaling, dist
            two_body_distance[n,two_body_counts[n]] = dist 
            two_body_scaling[n,two_body_counts[n]] = scaling
            two_body_counts[n] += 1

    return two_body_distance, two_body_scaling, two_body_counts


def generate_three_body(coordinates, nuclear_charge):

    n = 0

    indexes = dict()

    for i in range(1, max_id +1):
        for j in range(1, max_id +1):
            for k in range(j, max_id +1):

                print n, i, j, k

                key = (i, j, k)

                indexes[key] = n

                n += 1

    three_body_angle = np.zeros((n, THREE_BODY_MAX))
    three_body_scaling = np.zeros((n, THREE_BODY_MAX))
    three_body_counts = np.zeros(n, dtype = np.int32)


    print indexes
    print n

    for i, nuc1 in enumerate(nuclear_charge):

        q1 = nuc1 - 1
        for j, nuc2 in enumerate(nuclear_charge):

            if (i == j): continue
            q2 = nuc2 - 1

            for k, nuc3 in enumerate(nuclear_charge):

                if ((j <= k) or (k == i)): continue

                q3 = nuc2 - 1

                key = (nuc1, min(nuc2, nuc3), max(nuc2, nuc3))
                n = indexes[key]

                angle = bondangle(coordinates[j], coordinates[i], coordinates[k])
                
                angle2 = bondangle(coordinates[i], coordinates[j], coordinates[k])
                angle3 = bondangle(coordinates[i], coordinates[k], coordinates[j])

                dist1 = distance(coordinates[i], coordinates[j])
                dist2 = distance(coordinates[i], coordinates[k])
                dist3 = distance(coordinates[j], coordinates[k])

                scaling = (1 + 3.0 * np. cos(angle) * np.cos(angle2) * np.cos(angle2)) \
                            / (dist1*dist2*dist3)**2


                print i, j, k, key, n, angle / np.pi * 180.0, scaling


                three_body_angle[n,three_body_counts[n]] = angle
                three_body_scaling[n,three_body_counts[n]] = scaling
                three_body_counts[n] += 1

    print three_body_angle
    print three_body_scaling
    print three_body_counts


class BOSS(object):

    def __init__(self, coordinates, nuclear_charges):
        self.coordinates = coordinates
        self.nuclear_charges = nuclear_charges

        self.one_body, self.one_body_counts = generate_one_body(coordinates, nuclear_charges)

        (a, b, c) =  generate_two_body(coordinates, nuclear_charges)
        self.two_body_distance = a
        self.two_body_scaling = b
        self.two_body_counts = c

def distance(a,b):

    return np.sqrt(np.sum(np.square(a - b)))

def bondangle(a,b,c):
    # In case numpy.dot() returns larger than 1
    # and we cannot take acos() to that number

    acos_out_of_bound = 1.0
    v1 = a - b
    v2 = c - b
    v1 = v1 / np.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2)
    v2 = v2 / np.sqrt(v2[0]**2 + v2[1]**2 + v2[2]**2)
    dot_product = np.dot(v1,v2)

    if dot_product > acos_out_of_bound:
        dot_product = acos_out_of_bound
    if dot_product < -1.0 * acos_out_of_bound:
        dot_product = -1.0 * acos_out_of_bound

    return np.arccos(dot_product)


    rep = BOSS()


    self.coordinates = Nong
    self.one_body = None
    self.one_body_counts = None
    
    self.two_body_distance = None
    self.two_body_scaling = None
    self.two_body_counts = None
    return two_body_distances, two_body_scaling, two_body_counts

def kernel(boss1, boss2, sigma):

    # print boss1.two_body_counts
    # print boss2.two_body_counts

    K = 0.0 

    cut = 4*sigma

    inv_sigma = -1.0 / (sigma**2)

    for i, count1 in enumerate(boss1.two_body_counts):
        if count1 < 1: continue
        count2 = boss2.two_body_counts[i]

        if count2 < 1: continue

        for a in range(count1):
            for b in range(count2):

                # print boss1.two_body_distance[i,a], boss2.two_body_distance[i,b]
                mu = (boss1.two_body_distance[i,a] - boss2.two_body_distance[i,b])**2

                if mu > cut: continue
                # print boss1.two_body_scaling[i,a], boss2.two_body_scaling[i,b]
                scaling = boss1.two_body_scaling[i,a] * boss2.two_body_scaling[i,b]

                contrib = np.exp( mu * inv_sigma) * scaling
                K += contrib
                # print contrib, mu, scaling

        # print i, count1, count2

    return K

if __name__ == "__main__":

    test_dir = os.path.dirname(os.path.realpath(__file__))

    # Parse file containing PBE0/def2-TZVP heats of formation and xyz filenames
    data = get_energies(test_dir + "/hof_qm7.txt")

    # Generate a list of qml.Compound() objects
    mols = []

    np.set_printoptions(linewidth=66666)

    for xyz_file in sorted(data.keys()):

        # Initialize the qml.Compound() objects
        mol = qml.Compound(xyz=test_dir + "/qm7/" + xyz_file)


        # Associate a property (heat of formation) with the object
        mol.properties = data[xyz_file]

        # This is a Molecular Coulomb matrix sorted by row norm
        # mol.generate_coulomb_matrix()

        mols.append(mol)
   
        # rep = BOSS(mol.coordinates, mol.nuclear_charges)
        # boss.append(rep)

    # c = qml.Compound(xyz="water.xyz")
    # c2 = qml.Compound(xyz="water2.xyz")
    # rep1 = BOSS(c.coordinates, c.nuclear_charges)
    # rep2 = BOSS(c2.coordinates, c2.nuclear_charges)

    # Shuffle molecules
    np.random.seed(666)
    np.random.shuffle(mols)

    # Make training and test sets
    n_test  = 60
    n_train = 140

    training = mols[:n_train]
    test  = mols[-n_test:]

    sigma = 100.0
    llambda = 1e-10
    
    # List of properties
    Y = np.array([mol.properties for mol in training])
    Ys = np.array([mol.properties for mol in test])

    # List of representations
    mbtypes = get_slatm_mbtypes(np.array([mol.nuclear_charges for mol in mols]))
    for i, mol in enumerate(training): 
        mol.generate_slatm(mbtypes, local=False, rpower=6)
    for i, mol in enumerate(test): 
        mol.generate_slatm(mbtypes, local=False, rpower=6)

    X  = np.array([mol.representation for mol in training])
    Xs = np.array([mol.representation for mol in test])
    
    # Generate training Kernel
    K = gaussian_kernel(X, X, sigma)
    Ks = gaussian_kernel(X, Xs, sigma)
    
    K[np.diag_indices_from(K)] += llambda
    alpha = cho_solve(K,Y)

    # Calculate prediction kernel
    Yss = np.dot(Ks.transpose(), alpha)

    mae = np.mean(np.abs(Ys - Yss))

    print mae
    exit()

    # sigma = 0.1
    # llambda = 1e-9

    # X  = [BOSS(mol.coordinates, mol.nuclear_charges) for mol in training]
    # Xs = [BOSS(mol.coordinates, mol.nuclear_charges) for mol in test]
   
    # print X[0].two_body_counts
    # exit
    # K = np.zeros((len(X), len(X)))
    # for i, rep1 in enumerate(X): 
    #     for j, rep2 in enumerate(X):
    #         K[i,j] = kernel(rep1, rep2, sigma)
    #     print i, K[i,:]

    # Ks = np.zeros((len(X), len(Xs)))
    # for i, rep1 in enumerate(X): 
    #     for j, rep2 in enumerate(Xs):

    #         Ks[i,j] = kernel(rep1, rep2, sigma)
    #     print i, Ks[i,:]

    # K[np.diag_indices_from(K)] += llambda
    # alpha = cho_solve(K,Y)

    # # Calculate prediction kernel
    # Yss = np.dot(Ks.transpose(), alpha)

    # mae = np.mean(np.abs(Ys - Yss))

    # print mae
   
