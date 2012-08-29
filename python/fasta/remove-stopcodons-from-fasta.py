#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <fasta> ")
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['fasta'] = value
    
  if not args.has_key('fasta'):
    stderr( "fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('fasta') ):
    stderr( "fasta file does not exist." )
    show_help()
  
  return args


# =============================================================================
def parse_fasta(file):
  hash = {}
  fo = open(file)
  STOPCODONS = ["TAA", "TGA", "TAG"]
  id = ""
  for line in fo:
    line = line.strip()
    if line.startswith(">"):
      id = line[1:]
      if id.count(" ") > 0: id = id[:id.index(" ")]
      hash[id] = ""
    else:
      sequence = line.upper()
      i = 0
      while i < len(sequence):
        codon = sequence[i:i+3]
        if codon in STOPCODONS:
          hash[id] += "---"
        else:
          hash[id] += codon
        i += 3
  return hash

# =============================================================================
def replace_stop_codons(hash):
  for id, sequence in hash.iteritems():
    i = 0
    while i < len(sequence):
      codon = sequence[i:i+3]
      if codon in STOPCODONS:
        sequence[i:i+3] = "---"
      i += 3

  return hash


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  hash = parse_fasta(args['fasta'])
  width = 60
  for id, sequence in hash.iteritems():
    print ">" + id
    i = 0
    while i < len(sequence):
      part = sequence[i:min([len(sequence),i+60])]
      print part
      i += 60

# =============================================================================
args = handle_arguments()
main( args )

