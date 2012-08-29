#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> -i <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -b        blastout file (-m 8)" )
  stdout( " -i        file with the IDs to keep" )
  stdout( " " )
  
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hi:b:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()
  
  args = {}
  args['verbose'] = 0
  for key, value in keys:
    if key == '-b': args['in-blastout'] = value
    if key == '-i':  args['in-ids'] = value

  if not args.has_key('in-blastout'):
    stderr( "in-blastout file missing." )
    show_help()
  if not args.has_key('in-ids'):
    stderr( "in-ids file missing." )
    show_help()
    
  if not file_exists( args.get('in-blastout') ):
    stderr( "in-blastout file does not exist." )
    show_help()
  if not file_exists( args.get('in-ids') ):
    stderr( "in-ids file does not exist." )
    show_help()
  
  return args

# =============================================================================
def get_ids_to_remove( args ):
  """
  reads in the in-ids file and gathers all IDs to which
  the out fasta file will be reduced to.
  """
  fo = open( args.get('in-ids'), 'r' )
  ids = {}
  for line in fo:
    line = line.rstrip()
    ids[ line.replace('>','') ] = 1
  fo.close()
  return ids
  
  
# =============================================================================
def reduce_blastout( args, rmids ):
  """
  reads in in-fasta and creates out-fasta that only contains the records
  whose id is contained in the hash keepids.
  """
  
  retained = 0
  fo = open( args.get('in-blastout') )
  for line in fo:
    line = line.rstrip()
    if len(line) == 0: continue
    hid, qid = line.split("\t")[0:2]
    if rmids.has_key(hid) or rmids.has_key(qid): continue
    print line
    retained += 1
  fo.close()
    

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
rmids = get_ids_to_remove( args )
reduce_blastout( args, rmids )
