#!/usr/bin/python
import os, sys, string, anydbm
from low import *
from orthomcl import OrthoMCLCluster


# =============================================================================
def usage():
  print >> sys.stderr, "add significant BLAST hits (e.g. in-paralogs) to an existing orthomcl cluster.\n"
  print >> sys.stderr, "usage:      (1) " + sys.argv[0] + " noparalogs.orthomcl.out  blastout.add.dbm" 
  print >> sys.stderr, "         or (2) " + sys.argv[0] + " noparalogs.orthomcl.out  all.fasta  all.gg  all.blastout" 
  sys.exit(1)


def plausi():
  if len(sys.argv) != 3 and len(sys.argv) != 5: usage()
  return sys.argv[1:]


def read_gg(inGG):
  outHash, speciesArray = {}, []
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


def get_seq_lengths(file):
  lengthHash, id = {}, ""
  fo = open(file)
  for line in fo: 
    line = line.strip()
    if line.startswith(">"):
      id = line[1:]
      if id.count(" ") > 0: id = id[:id.index(" ")]
      lengthHash[id] = 0 
    else: lengthHash[id] += len(line)
  return lengthHash


def main():
  args = plausi()
  in_orthomcl = args[0]
  EVALUE = float('1e-20')
  IDENTITY = 30.0
  if len(args) == 4:
    in_fasta, in_gg, in_blast = args[1:4]
    gene2species, speciesArray = read_gg(in_gg)
    gene2length = get_seq_lengths(in_fasta)
    dbmfile = in_blast + ".add.dbm"
    dbm = anydbm.open(dbmfile, "c")
    fo = open(in_blast)
    for line in fo: 
      line = line.rstrip()
      cols = line.split("\t")
      qid, hid, evalue, identity = cols[0], cols[1], float(cols[10]), float(cols[2])
      # ignore self-hits and between-species hits, check e-value threshold
      if qid == hid: continue
      if gene2species[qid] != gene2species[hid]: continue
      if evalue > EVALUE: continue
      if identity < IDENTITY: continue
      # check that blast alignment spans at least 75% of the longer sequence
      alnlength, qlength, hlength = int(cols[3]), gene2length[qid], gene2length[hid]
      lengthcutoff = 0.80 * max([qlength, hlength])
      if alnlength < lengthcutoff: continue
      if not dbm.has_key(qid): dbm[qid] = ""
      else: dbm[qid] += " "
      dbm[qid] += hid
    fo.close()
    dbm.close()
  else: dbmfile = args[1]
  dbm = anydbm.open(dbmfile)

  fo = open(in_orthomcl)
  for line in fo:
    o = OrthoMCLCluster(line.rstrip())
    oldsize = o.get_count()
    additions = []
    for geneid, species in o.get_gene_hash().iteritems():
      if not dbm.has_key(geneid): continue
      [additions.append([x, species]) for x in dbm[geneid].split()]

    for x, species in additions: o.add_gene(x,species)
    o.to_s()
    newsize = o.get_count()
    print >> sys.stderr, "%s\t%s\t%s" %(o.get_name(), oldsize, newsize)

main()
