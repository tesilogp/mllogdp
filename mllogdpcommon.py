import numpy as np

def to_descriptor (orderedset, atomtypes):

    desc = np.zeros(len(orderedset))

    for at in atomtypes:
        idx = orderedset.index(at)

        desc[idx] += 1

    return desc

def import_descriptor (filename):
    fp = open(filename, "r")

    smilelogd = {}
    descriptors = []
    for l in fp:
        sline = l.split(",")
        smile = sline[0]
        logd = np.float64(sline[-1])

        descriptors.append([ np.int(v) for v in sline[1:-1]])

        smilelogd[smile] = logd
    
    fp.close()

    return smilelogd, np.asarray(descriptors)


