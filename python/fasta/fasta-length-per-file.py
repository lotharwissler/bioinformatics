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
  length = 0
  for line in fo:
    line = line.strip()
    if line.startswith(">"):
      continue
    else:
      length += len(line)
  fo.close()
  base = file
  if base.count(".") > 0: base = base[:base.index(".")]
  if base.count("_") > 0: base = base[:base.index("_")]
  print base + "\t" + str(length)


def main():
  inFile = plausi()
  parse_fasta_file( inFile )


main()
