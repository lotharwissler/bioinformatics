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
  stdout( " -f        fast file which should be reduced" )
  stdout( " -i        file with the IDs to keep" )
  stdout( " -v        verbose: report statistics to STDERR, otherwise silent" )
  stdout( " " )
  
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hi:f:v" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()
  
  args = {}
  args['verbose'] = 0
  for key, value in keys:
    if key == '-f': args['in-fasta'] = value
    if key == '-i':  args['in-ids'] = value
    if key == '-v':  args['verbose'] = 1

  if not args.has_key('in-fasta'):
    stderr( "in-fasta file missing." )
    show_help()
  if not args.has_key('in-ids'):
    stderr( "in-ids file missing." )
    show_help()
    
  if not file_exists( args.get('in-fasta') ):
    stderr( "in-fasta file does not exist." )
    show_help()
  if not file_exists( args.get('in-ids') ):
    stderr( "in-ids file does not exist." )
    show_help()
  
  return args

# =============================================================================
def get_ids_to_keep( args ):
  """
  reads in the in-ids file and gathers all IDs to which
  the out fasta file will be reduced to.
  """
  fo = open( args.get('in-ids'), 'r' )
  keepids = {}
  for line in fo:
    line = line.rstrip()
    keepids[ line.replace('>','') ] = 1
  fo.close()
  return keepids
  
  
# =============================================================================
def reduce_fasta( args, keepids ):
  """
  reads in in-fasta and creates out-fasta that only contains the records
  whose id is contained in the hash keepids.
  """
  if args.get('verbose'):
    sys.stderr.write('\tnumber of records to retain: %s ' % len(keepids) )
  retained = 0
  id, seq = "", ""
  fo = open( args.get('in-fasta') )
  for line in fo:
    line = line.rstrip()
    if len(line) == 0: continue
    if line[0] == ">":
      if id != "" and seq != "" and keepids.has_key(id):
        print ">" + id + "\n" + seq
        retained += 1
        id, seq = "", ""
      checkid = line[1:].split()[0]
      if keepids.has_key(checkid): id = checkid

    else:
      if id != "":
        if seq != "": seq += "\n"
        seq += line

  fo.close()
  if id != "" and seq != "" and keepids.has_key(id):
    print ">" + id + "\n" + seq
    retained += 1
  if args.get('vebose'):
    sys.stderr.write('| retained: %s | done.\n' % retained  )

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
keepids = get_ids_to_keep( args )
reduce_fasta( args, keepids )
