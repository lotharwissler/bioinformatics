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
  stdout( "usage: " + sys.argv[0] + " -f <path> -i -n" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        input file" )
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
    if key == '-f': args['file'] = value
    
  if not args.has_key('file'):
    stderr( "input file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "input file does not exist." )
    show_help()
  
  return args

  
# =============================================================================
def read_input(file):
  hash = {}
  speciesarray = []
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    pair, rate = line.split("\t")
    rate = str(round(1-float(rate),4))
    while len(rate) < 6: rate += "0"
    hash[pair] = rate
    speciesarray.extend(pair.split(","))
  fo.close()
  speciesarray = list(set(speciesarray))
  speciesarray.sort()
  return speciesarray, hash


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  speciesarray, hash = read_input(args['file'])
  print "\t" + str(len(speciesarray)+1)
  print "outgroup  " + "0.0000" + "\t" + string.join(["1.0000"]*len(speciesarray), "\t")
  for sp1 in speciesarray:
    line = sp1
    while len(line) < 10: line += " "
    line += "1.0000"
    for sp2 in speciesarray:
      key = [sp1,sp2]
      key.sort()
      key = string.join(key, ",")
      if sp1 == sp2: line += "\t" + "0.0000"
      else: line += "\t" + hash[key]
    print line

# =============================================================================
args = handle_arguments()
main( args )

