#!/usr/bin/python
import os, sys

def usage():
  print >> sys.stderr, "reports all fasta files with one or more sequences < or > n characters"
  print >> sys.stderr, "usage: " + sys.argv[0] + " folder \"<> n\""
  sys.exit(1)


def plausi():
  if len(sys.argv) != 3: usage()
  inFolder, inCutoff = sys.argv[1:3]
  inCut, inThreshold = inCutoff.split()
  inThreshold = int(inThreshold)
  return inFolder, inCut, inThreshold


def parse_fasta_file(file):
  lengthHash = {}
  fo = open(file)
  id = ""
  for line in fo:
    line = line.strip()
    if line.startswith(">"):
      id = line[1:]
      lengthHash[id] = 0
    else:
      lengthHash[id] += len(line)
  lengths = lengthHash.values()
  lengths.sort()
  return lengths[0]


def test_threshold(length, inCut, inThreshold):
  if inCut == ">" and length > inThreshold: return 1
  if inCut == "<" and length < inThreshold: return 1
  return 0


def main():
  inFolder, inCut, inThreshold = plausi()
  for filename in os.listdir(inFolder):
    if not filename.endswith(".fasta"): continue
    minlength = parse_fasta_file( filename )
    report = test_threshold(minlength, inCut, inThreshold)
    if report: 
      print os.path.split(filename)[1] + "\t" + str(minlength)
    


main()
