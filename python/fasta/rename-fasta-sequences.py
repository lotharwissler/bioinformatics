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
  stdout( " -f        fasta file to import" )
  stdout( " -m        tab delimited file that maps a regex to the replacement name, one per line" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:m:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-m': args['mapping'] = value
    
  if not args.has_key('file'):
    stderr( "fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "fasta file does not exist." )
    show_help()
    
  if not args.has_key('mapping'):
    stderr( "mapping file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "mapping file does not exist." )
    show_help()
 
  return args


# =============================================================================
def get_mapping(mfile):
  hash = {}
  fo = open( mfile, "r" )
  for line in fo:
    line = line.rstrip()
    if len(line) == 0: break
    if len(line.split("\t")) != 2: continue
    regex, replacement = line.split("\t")
    hash[re.compile(regex)] = replacement
  fo.close()
  return hash

# =============================================================================
def apply_replacement(idline, maphash):
  id = idline[1:].split()[0]
  for regex, replacement in maphash.iteritems():
    if re.search(regex, idline[1:]):
      idline = '>' + re.sub(regex, replacement, idline[1:], count=1)
      break
  return idline

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  maphash = get_mapping( args.get('mapping') )
  
  fo = open( args.get('file') )
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"): line = apply_replacement(line, maphash)
    print line
  fo.close()

# =============================================================================
args = handle_arguments()
main( args )

