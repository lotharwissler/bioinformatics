#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from goterm import GOTerm
from collections import defaultdict


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        DNA fasta file" )
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
    if key == '-f': args['file'] = value
    
  if not args.has_key('file'):
    stderr( "fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "fasta file does not exist." )
    show_help()
  
  return args

  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  counts = {'A':0, 'T':0, 'G':0, 'C':0}
  fo = open(args['file'])
  for line in fo:
    if line.startswith(">"): continue
    line = line.rstrip().upper()
    for char in ['A', 'T', 'G', 'C']:
      counts[char] += line.count(char)

  total = sum(counts.values())
  gc = 1.0 * (counts['G'] + counts['C']) / total
  base = args['file']
  if base.count(".") > 0: base = base[:base.index(".")]
  if base.count("_") > 0: base = base[:base.index("_")]

  print base + "\t" + str(gc)

# =============================================================================
args = handle_arguments()
main( args )

