#!/usr/bin/python
import os, sys, re

MODELS = ["M0", "M7", "M8", "Free", "M1a", "M2a", "MT1", "MT2", "MT3", "MT4", "MT5", "MT6"]


def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " folder"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 2: usage()
  inFolder = sys.argv[1]
  return inFolder


def get_all_base_files(inFolder):
  fileHash = {}
  for file in os.listdir(inFolder):
    filename = os.path.split(file)[1]
    basename = filename
    while basename.count('.') > 0: basename = os.path.splitext(basename)[0]
    fileHash[basename] = 1
  return fileHash.keys()


def parse_all_from_basefile(file):
  filesToParse = []
  for m in MODELS: filesToParse.append(file + ".codeml." + m)
  for f in filesToParse:
    if not os.path.exists(f) or not os.path.isfile(f): 
      print >> sys.stderr, "bad stuff happening with file", file, "/", f
      return

  modelHash = {}
  for f in filesToParse:
    fo = open(f)
    for line in fo:
      if line.startswith("lnL("):
        np = re.match("lnL.*\s+np:\s*(\d+)", line ).group(1)
        lnL = re.match("lnL\(.*\):\s+([0-9.-]+)", line ).group(1)
        break
    modelHash[ os.path.splitext(f)[1][1:] ] = [lnL, np]
    fo.close()

  sys.stdout.write(file)
  for m in MODELS:
    sys.stdout.write("\t" + m + ":" + modelHash[m][0] + "," + modelHash[m][1])
  sys.stdout.write("\n")


def main():
  inFolder = plausi()
  basefiles = get_all_base_files(inFolder)
  for basefile in basefiles:
    parse_all_from_basefile(basefile)


main()
