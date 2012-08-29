#!/usr/bin/python
import os, sys, re

def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " folder  (files end with *.MAalt)"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 2: usage()
  inFolder = sys.argv[1]
  return inFolder


def parse_from_file(inFile):
  basename = os.path.split(inFile)[1]
  fo = open(inFile)
  line = fo.readline().rstrip()
  while 1:
    if line.startswith("ns ="): 
      print >> sys.stderr, inFile
      length = re.search("ls =\s+(\d+)", line).group(1)

    if not line.startswith("Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)"):
      line = fo.readline().rstrip()
    else:
      line = fo.readline().rstrip() # Positive sites for foreground lineages Prob(w>1):
      line = fo.readline().rstrip()
      if re.match("^$", line): 
        sites = 0
      else:
        sites = 0
        while not re.match("^$", line):
          if line.endswith("*"):
            sites += 1
          print basename + "\t" + length + "\t" + line
          line = fo.readline().rstrip()
      break
  print >> sys.stderr, basename + "\t" + str(sites)
  fo.close()


def parse_all_files(inFolder):
  for filename in os.listdir(inFolder):
    if not filename.endswith(".MAalt"): continue
    parse_from_file(inFolder + "/" + filename)


def main():
  inFolder = plausi()
  parse_all_files(inFolder)


main()
