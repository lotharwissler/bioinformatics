#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import math        # match functions
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        GO flat file to import [tab delimited]" )
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
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
    

  return args


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  from collections import defaultdict
  goHash = defaultdict(list)
  fo = open( args.get('file') )
  for line in fo:
    line = line.strip()
    geneid, goterm = line.split("\t")
    if geneid.count(" ") > 0:
      geneid = geneid[:geneid.index(" ")]
    goHash[geneid].append(goterm)
  fo.close()
  for geneid, goterms in goHash.iteritems():
    print geneid + "\t" + string.join(goterms, ", ")

# =============================================================================
args = handle_arguments()
main( args )

