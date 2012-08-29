#!/usr/bin/python
import os, sys
from low import *

# takes an input protein fasta file and an orthomcl.gg file
# orthomcl.gg file format:
# speciesname1: id1 id2 id3 id4 .... full genome
# speciesname2: id1 id2 id3 id4 .... full genome
#
# with these infos, the goal is to get only one protein sequence per species
# we use t-coffee to find the most similar protein sequence per species
# to the whole cluster. so in case one species contributes several sequences 
# to a cluster, we choose the one species to keep which has the highest average 
# similarity to the rest of the cluster. if more than 1 sequence yield the highest
# avgsim, we determine whether these protein sequences are (1) all identical, 
# or whether they are (2) slightly different. In case (1), we choose any sequence
# randomly because it does not matter. In case (2), we sum up all pairwise
# similarities for each candidate sequence, and keep only the one sequence
# with the highest sum. If these are identical as well, we again choose randomly
# (should happen very rarely).



def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " clustering.out  orthomcl.gg"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 3: usage()
  inClustering, inGG = sys.argv[1:3]
  return inClustering, inGG


def get_number_of_species(inGG):
  count = 0
  fo = open(inGG)
  for line in fo: count += 1
  fo.close()
  return count


def read_gg(inGG):
  outHash = {}
  speciesArray = []
  fo = open(inGG)
  for line in fo:
    line = line.rstrip()
    cols = line.split()
    species = str(cols[0])[:-1]
    if not species in speciesArray: speciesArray.append(species)
    for col in cols[1:]:
      outHash[col] = species
  fo.close()
  return outHash, speciesArray


def main():
  inClustering, inGG = plausi()
  speciesHash, speciesArray = read_gg(inGG)
  
  fo = open(inClustering)
  for line in fo:
    if line.startswith("#"): continue
    line = line.rstrip()
    cluster, count, geneids = line.split("\t")[0:3]
    geneids = geneids.split(", ")
    currentSpecies = []
    for id in geneids: currentSpecies.append(speciesHash[id])
    speciesCount = len(set(currentSpecies))
    sys.stdout.write("%s(%s genes, %s taxa):\t" %(cluster, count, speciesCount)) 
    for id in geneids: 
      species = speciesHash[id]
      sys.stdout.write(id + "(" + species + ") ")
    sys.stdout.write("\n")
  fo.close()


main()
