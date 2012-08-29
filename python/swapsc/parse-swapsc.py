#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file" )
  stdout( " -b        branches of interest comma-separated" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:b:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f':	args['file'] = value
    if key == '-b':	args['branches'] = value.split(",")
        
  if not args.has_key('file'):
    stderr( "swapsc out file missing." )
    show_help()
  if not file_exists( args.get('file') ):
    stderr( "swapsc out file does not exist." )
    show_help()
    
  return args

# =============================================================================
def get_contraints(file):
  contraints = {}
  interest = 0
  fo = open(args['file'])
  for line in fo:
    if line.startswith("Proportion of codon sites under selective constraints"): break
    if not interest and not line.startswith("============================================================================================================="): continue
    if line.startswith("============================================================================================================="): 
      interest = 1
      continue
    line = line.rstrip()
    columns = line.split()
    if len(columns) == 1:
      currentbranch = columns[0]
      contraints[currentbranch] = {}
    else:
      if len(columns) < 9: continue
      if columns[7] == ">": continue # not significant signal
      contraints[currentbranch][columns[0]] = columns[9]
  for line in fo:
    line = line.rstrip()
    columns = line.split()
    if len(columns) != 5 or columns[0] != "S": continue
    S = columns[1]
  fo.close()
  return contraints, S

# =============================================================================
# =============================================================================
def main( args ):

  constraintshash, S = get_contraints(args['file'])
  outarray = [args['file'], S]
  for branch in args['branches']:
    if not constraintshash.has_key(branch): 
      outarray.append("0")
      continue
    PS = "0"
    for coord, type in constraintshash[branch].iteritems():
      if type == "PS": PS = "1"
    outarray.append(PS)
  print string.join(outarray, "\t")

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
