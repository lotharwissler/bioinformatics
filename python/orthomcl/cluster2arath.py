#!/usr/bin/python
import os, sys, string
from low import *
from orthomcl import OrthoMCLCluster


# =============================================================================
def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " noparalogs.orthomcl.out" 
  sys.exit(1)


def plausi():
  if len(sys.argv) != 2: usage()
  inFile = sys.argv[1]
  return inFile


def main():
  inFile = plausi()
  fo = open(inFile)
  for line in fo:
    o = OrthoMCLCluster(line.rstrip())
    print o.get_name() + "\t" + o.get_species_hash()['Arath'][0]


main()
