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
  stdout( " -f        file to map to" )
  stdout( " -m        file to map from" )
  stdout( " -a        column to index in the map-to file" )
  stdout( " -d        delimiter (default: tab | allowed: ; , tab space" )
  stdout( " -v        verbose/debug mode" )
  stdout( " -s        silent mode: no stderr msgs" )
  stdout( " -c        conservative: output only entries where a mapping was found" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hvscf:m:d:a:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  args['silent'] = 0
  for key, value in keys:
    if key == '-f': args['mapto'] = value
    if key == '-m': args['mapfrom'] = value
    if key == '-d': args['delimiter'] = value
    if key == '-a': args['colto'] = int(value)
    if key == '-v': args['debug'] = 1
    if key == '-s': args['silent'] = 1
    if key == '-c': args['conserve'] = 1
    
  if not args.has_key('mapto'):
    stderr( "map-to file argument missing." )
    show_help()
  elif not file_exists( args.get('mapto') ):
    stderr( "map-to file does not exist." )
    show_help()

  if not args.has_key('mapfrom'):
    stderr( "map-from file argument missing." )
    show_help()
  elif not file_exists( args.get('mapfrom') ):
    stderr( "map-from file does not exist." )
    show_help()
    
  if not args.has_key('delimiter') or args.get('delimiter') not in [ ";", ",", "tab", "space" ]: 
    args['delimiter'] = "\t"
  else:
    if args['delimiter'] == "tab": args['delimiter'] = "\t"
    elif args['delimiter'] == "space": args['delimiter'] = " "

  if not args.has_key('colto'): args['colto'] = 0
  if not args.has_key('debug'): args['debug'] = 0
  if not args.has_key('conserve'): args['conserve'] = 0

  return args


# =============================================================================
def get_mapping(file, delimiter, col, debug=0):
  hash = {}
  fo = open( file )
  for line in fo:
    line = line.rstrip()
    #print "delimiter:" + delimiter
    col = line.split(delimiter)
    if hash.has_key(col[0]): print >> sys.stderr, "*** WARNING: mapping ambiguous for entry", col[0]
    hash[ col[0] ] = string.join(col[1:], delimiter)
    if debug: print col[0], "=>", string.join(col[1:], delimiter)
  fo.close()
  return hash


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  mapping = get_mapping( args.get('mapfrom'), args.get('delimiter'), 0, args.get('debug'))
  fo = open( args.get('mapto') )
  for line in fo:
    line = line.rstrip()
    col = line.split(args.get('delimiter'))
    key = col[args.get('colto')]
    if mapping.has_key(key):
      col.append( mapping.get(key) )
    else:
      if not args.get('silent'): stderr('skipping entry. key not found in mapping-from file: \"%s\"' % key )
    if mapping.has_key(key) or not args.get('conserve'): 
      print string.join(col, args.get('delimiter'))

  fo.close()

# =============================================================================
args = handle_arguments()
main( args )

