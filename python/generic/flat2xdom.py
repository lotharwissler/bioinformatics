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
  stdout( " -p        prefix to put in fron of the key" )
  stdout( " -d        delimiter (default: space | allowed: ; , tab space" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:p:d:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-p': args['prefix'] = value
    if key == '-d': args['delimiter'] = value
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
    
  if not args.has_key('delimiter') or args.get('delimiter') not in [ ";", ",", "tab", "space" ]: 
    args['delimiter'] = 'space'

  return args


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  fo = open( args.get('file') )
  oldid = ""
  for line in fo:
    line = line.rstrip()
    if args.get('delimiter') == "tab":
      columns = line.split("\t")
    elif args.get('delimiter') == "space":
      columns = line.split()
    else:
      columns = line.split( args.get('delimiter') )
    id = columns[0]
    if id != oldid:
      oldid = id
      if args.has_key('prefix'):
        print ">" + args.get('prefix') + id
      else:
        print ">" + id
    print string.join( columns[1:], "\t" )
  fo.close()

# =============================================================================
args = handle_arguments()
main( args )

