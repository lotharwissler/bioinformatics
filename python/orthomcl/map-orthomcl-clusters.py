#!/usr/bin/python
import os, sys, string
from low import *


def usage():
  print >> sys.stderr, "usage: ", sys.argv[0], " <from>  <to>  [<orthomcl.out>]"
  print >> sys.stderr, "from/to: speciesname or \"cluster\""
  sys.exit(1)


def plausi():
  if len(sys.argv) < 3: usage()
  inTo = sys.argv[2].lower()
  inFrom = sys.argv[1].lower()
  if len(sys.argv) > 3:
    inTable = sys.argv[3]
  else:
    inTable = "/home/low/workspace/back-to-the-sea-orf-cluster-verification/32-new-orthologs/3-orthomcl-v1.4/noparalogs_orthomcl.out"
  if not os.path.exists(inTable) or not os.path.isfile(inTable) or not os.path.getsize(inTable) > 0: 
    print >> sys.stderr, "specified orthomcl table file does not exist, is not a file, or is empty\n"
    usage()
  return inFrom, inTo, inTable


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
  inFrom, inTo, inTable = plausi()
  fo = open(inTable)
  for line in fo:
    o = OrthoCluster(line.rstrip())
    speciesHash = o.get_species_hash()
    name = o.get_name()
    mapfrom, mapto = "", ""
    if inFrom == "cluster": mapfrom = name
    else: mapfrom = speciesHash[inFrom][0]
    if inTo == "cluster": mapto = name
    else: mapto = speciesHash[inTo][0]
    print mapfrom + "\t" + mapto
  fo.close()


main()
