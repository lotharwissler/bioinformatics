#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from collections import defaultdict
from misa import MisaSSRspecies


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        all.misa out file" )
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
  specieshash = {}
  fo = open(args['file'])
  for line in fo:
    m = MisaSSRspecies(line)
    if not specieshash.has_key(m.species): specieshash[m.species] = defaultdict(int)
    for char in ['A', 'T', 'G', 'C']:
      specieshash[m.species][char] += m.motif.count(char) * m.repeats

  speciesarray = specieshash.keys()
  speciesarray.sort()
  for species in speciesarray:
    total = sum(specieshash[species].values())
    gc = 1.0 * (specieshash[species]['G'] + specieshash[species]['C']) / total
    print species + "\t" + str(gc)

# =============================================================================
args = handle_arguments()
main( args )

