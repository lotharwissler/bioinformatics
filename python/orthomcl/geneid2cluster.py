#!/usr/bin/python
import os, sys, string
from low import *
from orthomcl import OrthoMCLCluster


# =============================================================================
def usage():
  print >> sys.stderr, "prints a mapping between each gene id and its cluster from orthomcl output\n"
  print >> sys.stderr, "usage: " + sys.argv[0] + " orthomcl.out" 
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
    name = o.get_name()
    geneHash = o.get_gene_hash()
    for geneid, species in geneHash.iteritems(): print geneid + "\t" + name


main()
