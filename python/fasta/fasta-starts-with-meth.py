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
  fo = open(file)
  id = ""
  for line in fo:
    line = line.strip()
    if line.startswith(">"):
      id = line[1:]
      if id.count(" ") > 0: id = id[:id.index(" ")]
      seq = ''
    elif len(seq) == 0:
      seq += line
      if seq[0].upper() == 'M': print "%s\t1" % id
      else: print "%s\t0" % id


def main():
  inFile = plausi()
  parse_fasta_file( inFile )


main()
