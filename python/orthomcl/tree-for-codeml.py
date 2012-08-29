#!/usr/bin/python
import os, sys
import string


def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " orthomcl.out  base.tree"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 3: usage()
  inOrtho, inTree = sys.argv[1:3]
  return inOrtho, inTree


class OrthoCluster():
  def __init__(self, line):
    descr, genedefs = line.split("\t")
    genedefs = genedefs.split()
    self.name = descr[:descr.index('(')].lower()
    self.geneHash = {}
    self.speciesHash = {}
    for genedef in genedefs:
      geneid = genedef[:genedef.index('(')]
      species = genedef[genedef.index('(')+1:-1] + "1"
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


def parse_orthocml_out(inFile, tree):
  fo = open(inFile)
  for line in fo:
    o = OrthoCluster(line.rstrip())
    speciesHash = o.get_species_hash()
    name = o.get_name()
    for species, genelist in speciesHash.iteritems():
      if len(genelist) > 1: break

    replacement = '(' + species[:-1] + '1 #1,' + species[:-1] + '2)'
    tree_repl_1 = tree.replace(species, replacement)
    replacement = '(' + species[:-1] + '1,' + species[:-1] + '2 #1)'
    tree_repl_2 = tree.replace(species, replacement)
    fw = open(name + ".tree.1", "w")
    fw.write(tree_repl_1)
    fw.close()
    fw = open(name + ".tree.2", "w")
    fw.write(tree_repl_2)
    fw.close()
  fo.close()


def read_tree_from_file(file):
  fo = open(file)
  tree = ""
  for line in fo:
    tree += line.strip()
  fo.close()
  return tree


def main():
  inOrtho, inTree = plausi()
  tree = read_tree_from_file(inTree)
  parse_orthocml_out(inOrtho, tree)



main()
