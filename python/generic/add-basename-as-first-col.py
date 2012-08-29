#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from collections import defaultdict
import fileinput


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        input file (will be rewritten on the fly!) - basename is everything before the first dot" )
  stdout( " -l        basename to lower case" )
  stdout( " -u        basename to upper case" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:ul" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'lower':False, 'upper':False}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-l': args['lower'] = True
    if key == '-u': args['upper'] = True
    
  if not args.has_key('file'):
    stderr( "fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "fasta file does not exist." )
    show_help()

  if args['lower'] and args['upper']:
    stderr( "cannot select both lower and upper." )
    show_help()
  
  return args

  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  filename = os.path.split(args['file'])[1]
  basename = filename
  while basename.count(".") > 0: basename = os.path.splitext(basename)[0]
  if args['lower']: basename = basename.lower()
  if args['upper']: basename = basename.upper()
  for line in fileinput.input(args['file'],inplace=1):
    print basename + "\t" + line.rstrip()

# =============================================================================
args = handle_arguments()
main( args )

