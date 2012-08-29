#!/usr/bin/python
import os, sys, string
from low import *


def usage():
  print >> sys.stderr, "usage: ", sys.argv[0], " [<noparalogs.orthomcl.out>]"
  sys.exit(1)


def plausi():
  if len(sys.argv) < 2: usage()
  inOrtho = sys.argv[1]
  if not os.path.exists(inOrtho) or not os.path.isfile(inOrtho) or not os.path.getsize(inOrtho) > 0: 
    print >> sys.stderr, "specified orthomcl file does not exist, is not a file, or is empty\n"
    usage()
  return inOrtho


class OrthoCluster():
  def __init__(self, line):
    descr, genedefs = line.split("\t")
    genedefs = genedefs.split()
    self.name = descr[:descr.index('(')].lower()
    self.geneHash = {}
    self.speciesHash = {}
    for genedef in genedefs:
      geneid = genedef[:genedef.index('(')]
      species = genedef[genedef.index('(')+1:-1].lower()
      self.geneHash[geneid] = species
      if self.speciesHash.has_key(species): self.speciesHash[species].append(geneid)
      else: self.speciesHash[species] = [geneid]

  def get_name(self): return self.name
  def get_count(self): return len(self.geneHash)
  def get_gene_hash(self): return self.geneHash
  def get_species_hash(self): return self.speciesHash
    


def main():
  inOrtho = plausi()
  fo = open(inOrtho)
  speciesCols = 0
  for line in fo:
    o = OrthoCluster(line.rstrip())
    SH = o.get_species_hash()
    if not speciesCols:
      speciesCols = SH.keys()
      speciesCols.sort()
      print "OrthoMCL.ID" + "\t" + string.join(speciesCols, "\t")

    name = o.get_name()
    print name + "\t" + string.join( [SH[x][0] for x in speciesCols], "\t")

  fo.close()


main()
