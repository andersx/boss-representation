# MIT License (c) Anders Christensen 2018

import numpy as np

def generate_input(nuclear_charges, coordinates, max_size=23):

    rep = np.zeros((4,max_size))
    rep[0,:] = -1.0

    for i, data in enumerate(zip(nuclear_charges, coordinates)):

        rep[0,i] = float(data[0])
        rep[1:4,i] = data[1][:3]

    return rep
