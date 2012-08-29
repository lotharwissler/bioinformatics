#!/usr/bin/python
import os, sys, string
import anydbm         # index databases (file hash)
from low import *
from collections import defaultdict

FASTAID_REGEX = re.compile(">(\S+)")

def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + "  orthomcl.out  all.fasta"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 3: usage()
  inFasta = sys.argv[2]
  inOrtho = sys.argv[1]
  if not os.path.exists(inFasta) or not os.path.isfile(inFasta) or not os.path.getsize(inFasta) > 0: 
    print >> sys.stderr, "specified input fasta file does not exist, is not a file, or is empty\n"
    usage()
  if not os.path.exists(inOrtho) or not os.path.isfile(inOrtho) or not os.path.getsize(inOrtho) > 0: 
    print >> sys.stderr, "specified input orthomcl.out file does not exist, is not a file, or is empty\n"
    usage()
  return inOrtho, inFasta


def cache_genomes(file, recreate=0):
  outdbm = file + ".dbm"
  if os.path.exists(outdbm) and os.path.getsize(outdbm) > 0 and not recreate: return outdbm
  DBM = anydbm.open( outdbm, 'c' )
  fo = open(file)
  key = ""
  for line in fo:
    line = line.strip()
    if line.startswith(">"): 
      key = re.match(FASTAID_REGEX, line).group(1)
      DBM[key] = ""
    else:
      DBM[key] += line
  DBM.close()
  return outdbm


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
      if species == "Osa":
        geneid = geneid[:geneid.index('|protein')] + '|CDS'
      elif species == "Vitis":
        geneid = "GSVIVT" + geneid[6:]
      elif species == "Zma":
        geneid = geneid.replace("_FGP", "_FGT")
        geneid = geneid.replace("_P", "_T")
      self.geneHash[geneid] = species
      if self.speciesHash.has_key(species): self.speciesHash[species].append(geneid)
      else: self.speciesHash[species] = [geneid]

  def get_name(self): return self.name
  def get_count(self): return len(self.geneHash)
  def get_gene_hash(self): return self.geneHash
  def get_species_hash(self): return self.speciesHash
    


def main():
  inOrtho, inFasta = plausi()
  info("    caching genomes ...")
  dbm = cache_genomes(inFasta)
  info(" done.")
  sout, serr = catch_bash_cmd_output( "wc -l %s" % inOrtho )
  total = int( sout.split()[0] )
  count = 0
  fo = open(inOrtho)
  for line in fo:
    o = OrthoCluster(line.rstrip())
    geneHash = o.get_gene_hash()
    name = o.get_name()
    idfile = name + ".ids"
    fastafile = name + ".fasta"
    ufastafile = name + ".ufasta"
    fw = open(idfile, 'w')
    for id, species in geneHash.iteritems(): fw.write(id + "\n")
    fw.close()
    
    fu = open(ufastafile, 'w')
    fw = open(fastafile, 'w')
    seqHash = anydbm.open(dbm, "r")
    collectedSpecies = defaultdict(int)
    for geneid, species in geneHash.iteritems():
      if not seqHash.has_key(geneid):
        print "ID", geneid, "not found"
        continue
      sequence = seqHash[geneid]
      collectedSpecies[species] += 1
      count = collectedSpecies[species]
      fw.write(">" + geneid + "\n")
      fu.write(">" + species + str(count) + "\n")
      i = 0
      width = 60
      while i < len(sequence):
        frac = sequence[i:min([len(sequence),i+width])]
        fw.write(frac + "\n")
        fu.write(frac + "\n")
        i += width

    seqHash.close()
    fw.close()
    fu.close()
    count += 1
    progress = int(50.0 * count / total) * "#"
    progress += (50 - len(progress)) * " "
    info("       0% " + progress + " 100%     ")

  info("       0% " + progress + " 100%    \n")
    


main()
