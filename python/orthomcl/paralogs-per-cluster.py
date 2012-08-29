#!/usr/bin/python
import os, sys
import string


def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " orthomcl.out"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 2: usage()
  inFile = sys.argv[1]
  return inFile


class OrthoCluster():
  def __init__(self, line):
    descr, genedefs = line.split("\t")
    genedefs = genedefs.split()
    self.name = descr[:descr.index('(')].lower()
    self.geneHash = {}
    self.speciesHash = {}
    for genedef in genedefs:
      geneid = genedef[:genedef.index('(')]
      species = genedef[genedef.index('(')+1:-1]
      self.geneHash[geneid] = species
      if self.speciesHash.has_key(species): self.speciesHash[species].append(geneid)
      else: self.speciesHash[species] = [geneid]

  def get_name(self): return self.name
  def get_count(self): return len(self.geneHash)
  def get_gene_hash(self): return self.geneHash
  def get_species_hash(self): return self.speciesHash


def get_species_from_first_line(inFile):
  fo = open(inFile)
  line = fo.readline()
  o = OrthoCluster(line.rstrip())
  fo.close()
  species = o.get_species_hash().keys()
  species.sort()
  return species


def parse_orthocml_out(inFile):
  speciesList = get_species_from_first_line(inFile)
  print >> sys.stdout, "\t" + string.join(speciesList, "\t")
  fo = open(inFile)
  for line in fo:
    o = OrthoCluster(line.rstrip())
    speciesHash = o.get_species_hash()
    sys.stdout.write(o.get_name())
    for s in speciesList:
      count = 0
      if speciesHash.has_key(s): count = len(speciesHash[s])
      sys.stdout.write("\t%s" % count)
    sys.stdout.write("\n")
        
  fo.close()


def main():
  inFile = plausi()
  parse_orthocml_out(inFile)



main()
