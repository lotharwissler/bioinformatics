#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from misa import MisaSSR
from collections import defaultdict


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        misa outptu file" )
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
    stderr( "misa file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "misa file does not exist." )
    show_help()
  
  return args


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  fo = open(args['file'])
  for line in fo:
    if line.startswith("ID\t"): continue
    m = MisaSSR(line)
    if m.type != "c" and m.type != "c*": print m.to_s()
    else:
      startpos = m.startpos
      separatepatterns = re.findall("\([ATGC]+\)\d+[*]{0,1}",m.pattern)
      for separatepattern in separatepatterns:
        motif = separatepattern[1:separatepattern.index(")")]
        if separatepattern.endswith("*"): repeats = int(separatepattern[separatepattern.index(")")+1:-1])
        else: repeats = int(separatepattern[separatepattern.index(")")+1:])
        length = len(motif)*repeats
        endpos = startpos + length -1
        print string.join([m.geneid, str(m.ssrnr), "p" + str(len(motif)), separatepattern, str(length), str(startpos), str(endpos)], "\t")
        startpos = endpos+1



# =============================================================================
args = handle_arguments()
main( args )

