#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself
import anydbm

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        nt alignment file" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:t:p:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f':	args['aln'] = value
        
  if not args.has_key('aln'):
    stderr( "aln file missing." )
    show_help()
  if not file_exists( args.get('aln') ):
    stderr( "aln file does not exist." )
    show_help()
    
  return args

# =============================================================================
def get_aln_length_from_file( filename ):
  fo = open( filename )
  firstline = fo.readline()
  n, length = firstline.split()
  fo.close()
  return length

# =============================================================================
def get_lnL_from_file( filename, model ):
  file = filename + '.paml.out.' + model
  np, lnL = None, None
  if not file_exists( file ): 
    stderr( "File does not exist: %s" %file )
    return np, lnL

  fo = open( file )
  for line in fo:
    if line.startswith("lnL"):
      #print filename, model, line
      np = re.match("lnL\(.*\s+np:\s*(\d+)", line ).group(1)
      lnL = re.match("lnL\(.*\):\s+([0-9.-]+)", line ).group(1)
      break
  fo.close()
  return np, lnL

# =============================================================================
# =============================================================================
def main( args ):
  
  models = ["M0", "M3K2", "M3K3", "M7", "M8", "Free"]
  filename = args.get('aln')
  
  line = []
  line.append( filename )
  length = get_aln_length_from_file( filename )
  line.append( length )
  for M in models:
    np, lnL = get_lnL_from_file( filename, M )
    if np == None or lnL == None:
      stderr( "%s: None returned for model %s (%s/%s)" %( filename, M, np, lnL ) )
      sys.exit(1)
    line.append( M )
    line.append( np )
    line.append( lnL )
  print string.join(line,"\t")

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
