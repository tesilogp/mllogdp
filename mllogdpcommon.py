import numpy as np

def to_descriptor (orderedset, atomtypes):

    desc = np.zeros(len(orderedset))

    for at in atomtypes:
        idx = orderedset.index(at)

        desc[idx] += 1

    return desc


