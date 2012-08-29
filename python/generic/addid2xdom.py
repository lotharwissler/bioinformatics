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
  stdout( "usage: " + sys.argv[0] + " -f <path> -i -n" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file to import" )
  stdout( " -i        id mapping file" )
  stdout( " -n        column to look up the id for [0..n]" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:i:n:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-i': args['idfile'] = value
    if key == '-n': args['column'] = int(value)
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
  
  if not args.has_key('idfile'):
    stderr( "import id file argument missing." )
    show_help()
  elif not file_exists( args.get('idfile') ):
    stderr( "import id file does not exist." )
    show_help()
  
  if not args.has_key('column'):
    stderr( "column argument missing." )
    show_help()

  return args


  return idhash
  
# =============================================================================
def get_idhash( args ):
  idhash = {}
  fo = open( args.get('idfile') )
  for line in fo:
    line = line.rstrip()
    if len(line.split("\t")) > 2:
      key = line.split("\t")[0]
      value = string.join(line.split("\t")[1:],"\t")
    else:
      key, value = line.split("\t")
    idhash[ key ] = value
  fo.close()
  return idhash

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  idhash = get_idhash( args )

  fo = open( args.get('file') )
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"): 
      id = line[1:]
      print line
      #asdf = 0
    else:
      columns = line.split("\t")
      if args.get('column') == 0:
        lookup = id
      else: 
        lookup = columns[ args.get('column') -1 ]
      if not idhash.has_key( lookup ):
        stderr( "lookup name not found in the id file: " + lookup )
        continue
        #sys.exit(1)
      id = idhash.get( lookup )
      columns.append(id)
      print string.join( columns, "\t" )
  fo.close()

# =============================================================================
args = handle_arguments()
main( args )

