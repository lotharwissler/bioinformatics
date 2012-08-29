#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from goterm import GOTerm
from collections import defaultdict


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        go obo file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['obo'] = value
    
  if not args.has_key('obo'):
    stderr( "obo file argument missing." )
    show_help()
  elif not file_exists( args.get('obo') ):
    stderr( "obo file does not exist." )
    show_help()
  
  return args

  
# =============================================================================
def read_obo( file ):
  hash = {}
  goterm = {}
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    if line.startswith("[Term]") or line.startswith("[Typedef]"):
      if goterm.has_key('id') and goterm.has_key('name'): hash[goterm['id']] = goterm['name']
      goterm = {}
    elif line.startswith("id:"):
      goterm['id'] = line.split()[1]
    elif line.startswith("name:"):
      goterm['name'] = string.join(line.split()[1:], " ")
  fo.close()
  print >> sys.stderr, "goterms read from obo: %s" % len(hash)
  return hash


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  gohash = read_obo(args['obo'])
  for goid, goname in gohash.iteritems():
    print goid + "\t" + goname

# =============================================================================
args = handle_arguments()
main( args )

