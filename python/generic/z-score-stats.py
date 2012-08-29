#!/usr/bin/python
import os, sys, string
from low import *
from collections import defaultdict
import rpy2.robjects as robjects
R = robjects.r


# =============================================================================
def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " paralog-count.tab" 
  sys.exit(1)


def plausi():
  if len(sys.argv) != 2: usage()
  inCounts = sys.argv[1]
  return inCounts


def R_mean_and_sd(pylist):
  rcountsvec = robjects.IntVector(pylist)
  mean = R['mean'](rcountsvec)[0]
  sd = R['sd'](rcountsvec)[0]
  return mean, sd


def Zscore(x, mean, sd):
  if sd == 0: return 0
  return (1.0*x - mean)/sd

def main():
  inCounts = plausi()
  fo = open(inCounts)
  lines = fo.readlines()
  fo.close()
  header = lines.pop(0).rstrip().split("\t")
  speciesArray = header[1:]
  results = defaultdict(lambda: defaultdict(int))
  for line in lines:
    line = line.rstrip()
    columns = line.split("\t")
    cluster = columns[0]
    genecounts = columns[1:]
    mean, sd = R_mean_and_sd(genecounts)
    for i in range(len(genecounts)):
      gc, species = int(genecounts[i]), speciesArray[i]
      z = Zscore(gc, mean, sd)
      if abs(z) < 2: continue
      if z > 3: results[species]['Z > 3'] += 1
      elif z > 2: results[species]['Z > 2'] += 1
      elif z < -3: results[species]['Z < -3'] += 1
      elif z < -2: results[species]['Z < -2'] += 1
  
  speciesArray.sort()
  print "\t" + string.join(speciesArray, "\t")
  for zcat in ['Z > 3', 'Z > 2', 'Z < -3', 'Z < -2']:
    sys.stdout.write(zcat)
    for spec in speciesArray:
      count = str(results[spec][zcat])
      sys.stdout.write("\t" + count)
    sys.stdout.write("\n")


main()
