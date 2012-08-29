#!/usr/bin/python

import sets
import sys, os

def get_lines_in_hash(file):
  hash = {}
  fo = open(file)
  for line in fo: hash[line.strip()] = 1
  fo.close()
  return hash

def get_lines( file ):
  lines = []
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    lines.append(line)

  return sets.Set(lines)

def terminate():
  print >> sys.stderr, "provide at least two valid input files as input arguments"
  sys.exit(1)


if len(sys.argv[1:]) < 2: terminate()
for inputfile in sys.argv[1:]:
  if not os.path.isfile(inputfile): terminate()

allhashes = []
for file in sys.argv[1:]:
  allhashes.append( get_lines_in_hash(file) )

refkeys = allhashes[0].keys()
for refkey in refkeys:
  found = 0
  for hash in allhashes:
    if hash.has_key(refkey): found += 1
    else: break
  if found == len(allhashes):
    print refkey

#l1 = get_lines(sys.argv[1])
#l2 = get_lines(sys.argv[2])
#for e in l1.intersection(l2):
#  print e
