#!/usr/bin/python
import os, sys

def usage():
  print >> sys.stderr, "reports for each fasta sequence the length in tab format"
  print >> sys.stderr, "usage: " + sys.argv[0] + " fastafile"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 2: usage()
  inFile = sys.argv[1]
  return inFile


def parse_fasta_file(file):
  lengthHash = {}
  fo = open(file)
  id = ""
  for line in fo:
    line = line.strip()
    if line.startswith(">"):
      id = line[1:]
      if id.count(" ") > 0: id = id[:id.index(" ")]
      lengthHash[id] = 0
    else:
      lengthHash[id] += len(line)
  for id, length in lengthHash.iteritems():
    print id + "\t" + str(length)


def main():
  inFile = plausi()
  parse_fasta_file( inFile )


main()
