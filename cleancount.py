# -*- coding: cp1252 -*-
import numpy as np

######################################

def get_homomsets(filename):
  homomsets=[]
  with open(filename,"r") as f:
      for i in range(8):
          f.readline()
      l = f.readline()
      while(not("found" in l)):
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

p=7

for n in range(8,16):
  homomsets_left=get_homomsets("./Card-"+str(p)+"/left-n"+str(n)+"p"+str(p)+".txt")
  homomsets_right=get_homomsets("./Card-"+str(p)+"/right-n"+str(n)+"p"+str(p)+".txt")


  print "N={} - Cardinality={}".format(n,p)
  val,counts = np.unique([len(x) for x in homomsets_left],return_counts=True)
  print "Left: ",zip(val,counts)
  val,counts = np.unique([len(x) for x in homomsets_right],return_counts=True)
  print "Right: ",zip(val,counts)

  c=[]
  for x in homomsets_left:
    if set(x) in [set(y) for y in homomsets_right]:
      c.append(len(x))
  val,counts = np.unique(c,return_counts=True)
  print "Left and Right homometric: ",zip(val,counts)

