import sys
import numpy as np

h=[]
hboth=[]
with open(sys.argv[1],"r") as f:
    t=f.readline()
    while(not t==""):
        t=f.readline()
        hashs = t.rstrip().replace(" ","").split("-")
        sets = f.readline(),f.readline()
        sets = list(zip(hashs,sets))
        t = f.readline()
        flag=0
        for y in h:
            if (sets[0] in y):
                if not sets[1] in y:
                    y.append(sets[1])
                flag=1
                break
            elif (sets[1] in y):
                if not sets[0] in y:
                    y.append(sets[0])
                flag=1
                break
        if not flag:
            h.append(sets)


        if "Right homometric" in t:
            flag=0
            for y in hboth:
                if (sets[0] in y):
                    if not sets[1] in y:
                        y.append(sets[1])
                    flag=1
                    break
                elif (sets[1] in y):
                    if not sets[0] in y:
                        y.append(sets[0])
                    flag=1
                    break
            if not flag:
                hboth.append(sets)
        t=f.readline()

print("# of left homometric subsets")
vals,counts = np.unique([len(x) for x in h],return_counts=True)
for x,y in zip(vals,counts):
    print("{}-uples: {}".format(x,y))
print("# of simultaneous left and right homometric subsets")
vals,counts = np.unique([len(x) for x in hboth],return_counts=True)
for x,y in zip(vals,counts):
    print("{}-uples: {}".format(x,y))
