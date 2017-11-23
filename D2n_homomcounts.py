# -*- coding: cp1252 -*-
import numpy as np

######################################

def get_homomsets(filename):
    """Quick an dirty function for reading an output file from the C
    enumeration program and listing the homometric sets.

    Parameters
    ----------
    filename : the name of the output file from the C program

    Returns
    -------
    The list of unique homometric n-uples
    """
    homomsets=[]
    with open(filename,"r") as f:
        l = f.readline()
        while(l):
            x,y = l.rstrip().split("-")
            x=int(x)
            y=int(y)
            f.readline()
            f.readline()
            f.readline()
            l = f.readline()
            flag=0
            for i,theset in enumerate(homomsets):
                if (x in theset) and not (y in theset):
                    flag=1
                    homomsets[i].append(y)
                if (y in theset) and not (x in theset):
                    flag=1
                    homomsets[i].append(x)
                if (y in theset) and (x in theset):
                  flag=1
            if not flag:
                homomsets.append([x,y])
    return homomsets


######################################

homomsets = get_homomsets("./output.txt")
print "Found {} unique homometric sets:".format(len(homomsets))
print homomsets
print "Counts by n-uples:"
val,counts = np.unique([len(x) for x in homomsets],return_counts=True)
for a,b in zip(val,counts):
    print "  Number of {}-uples: {}".format(a,b)
