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
  stdout( " -f        nt alignment file (fasta)" )
  stdout( " -m        paml M0 out file" )
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
    if key == '-f':	args['aln'] = value
    if key == '-m':	args['m0'] = value
        
  if not args.has_key('aln'):
    stderr( "aln file missing." )
    show_help()
  if not file_exists( args.get('aln') ):
    stderr( "aln file does not exist." )
    show_help()
    
  if not args.has_key('m0'):
    stderr( "M0 file missing." )
    show_help()
  if not file_exists( args.get('m0') ):
    stderr( "M0 file does not exist." )
    show_help()

  return args

# =============================================================================
# =============================================================================
def main( args ):

  #sys.stderr.write(args.get('aln') + "\t")
  #sys.stderr.flush()
  # create evolver control file based on the M0 out file
  fo = open( args.get('m0') )
  line = ""
  while not re.match("\s+\d+\s+\d+\s*$", line):
    line = fo.readline()
  numbers = line.split()
  nspecies, length = numbers[0:2] 
  fo.close()

  fo = open( args.get('aln') )
  print "  " + nspecies + "  " + length + "\n"
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"): print line[1:]
    else: print line
  fo.close()

  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
