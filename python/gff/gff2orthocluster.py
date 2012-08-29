#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
import gff3
from collections import defaultdict


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        gff3 file" )
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
    if key == '-f': args['gff'] = value
    
  if not args.has_key('gff'):
    stderr( "gff argument missing." )
    show_help()
  elif not file_exists( args.get('gff') ):
    stderr( "gff does not exist." )
    show_help()
  
  return args

  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  fo = open(args['gff'])
  for line in fo:
    if line.startswith("#"): continue
    if len(line.strip()) == 0: continue
    if len(line.split("\t")) != 9: continue
    gf = gff3.GeneFeature(line.rstrip())
    if gf.type != "gene": continue
    id = gf.get_attributes()['ID']
    if gf.strand == '+': strand = '1'
    else: strand = "-1"
    print string.join([id, gf.seqid, str(gf.start), str(gf.stop), strand], "\t")
  fo.close()

# =============================================================================
args = handle_arguments()
main( args )

