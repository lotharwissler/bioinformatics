#!/usr/bin/python
import os, sys, string


def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " folder with genomes (*.fasta or *.fasta.gz)"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 2: usage()
  inFolder = sys.argv[1]
  if not os.path.exists(inFolder) or not os.path.isdir(inFolder): 
    print >> sys.stderr, "specified input folder does not exist or is not a directory\n"
    usage()
  if not inFolder.endswith('/'): inFolder += '/'
  return inFolder


def iterate_folder(inFolder):
  inFiles = []
  for fname in os.listdir(inFolder):
    if not fname.endswith('.fasta') and not fname.endswith('.fasta.gz'): continue
    inFiles.append(inFolder + fname)
  return inFiles


def process_file(inFile):
  gzip = 0
  if inFile.endswith('.gz'): gzip = 1

  if gzip:
    ec = os.system('gunzip ' + inFile)
    inFile = os.path.splitext(inFile)[0]

  filename = os.path.split(inFile)[1]
  outName = os.path.splitext(filename)[0]

  sys.stdout.write(outName + ": ")

  ids = {}
  fo = open(inFile)
  for line in fo:
    if not line.startswith(">"): continue
    line = line.rstrip()
    id = line[1:]
    if id.count(" ") > 0: id = id[:id.index(" ")]
    ids[id] = 1

  sys.stdout.write( string.join(ids.keys(), " ") )
  sys.stdout.write("\n")

  if gzip: ec = os.system('gzip ' + inFile)
  

def main():
  inFolder = plausi()
  inFiles = iterate_folder(inFolder)
  for inFile in inFiles: process_file(inFile)

main()
