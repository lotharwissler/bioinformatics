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
  stdout( " -f        flat file to import [tab delimited]" )
  stdout( " -a        index of the first dimension key [default: 0]" )
  stdout( " -b        index of the second dimension key [default: 1]" )
  stdout( " -v        index of the value [default: 2]" )
  stdout( " -o        order: comma-separated list of keys in which to output the matrix [default: alphabetically sorted]" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:a:b:v:o:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'key1':0, 'key2':1, 'value':2}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-a': args['key1'] = int(value)
    if key == '-b': args['key2'] = int(value)
    if key == '-v': args['value'] = int(value)
    if key == '-o': args['order'] = value.split(",")
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
    

  return args


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  hash = {}
  keys = []
  fo = open( args.get('file') )
  for line in fo:
    col = line.strip().split("\t")
    key1, key2, value = col[args['key1']], col[args['key2']], col[args['value']]
    hash[key1 + '|||' + key2] = value
    if not key1 in keys: keys.append(key1)
    if not key2 in keys: keys.append(key2)
  fo.close()
  if args.has_key('order'): keys = args['order']
  else: keys.sort()

  print string.join(keys, ",")
  for i in keys:
    sys.stdout.write(i)
    for j in keys:
      value = 'NA'
      if hash.has_key(i+'|||'+j): value = hash[i+'|||'+j]
      elif hash.has_key(j+'|||'+i): value = hash[j+'|||'+i]
      sys.stdout.write(","+value)
    sys.stdout.write("\n")


# =============================================================================
args = handle_arguments()
main( args )

