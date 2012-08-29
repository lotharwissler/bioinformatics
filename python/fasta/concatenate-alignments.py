#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from collections import defaultdict


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> -i -n" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -e        file extension, e.g. \".muscle\"" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "he:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-e': args['ext'] = value
    
  if not args.has_key('ext'):
    stderr( "ext argument missing." )
    show_help()
  
  return args

  
# =============================================================================
def aln_is_conserved(file, min=0.85):
  popenout = os.popen("~/bin/t-coffee -other_pg seq_reformat -in %s -output sim | tail -n 1" % file)
  out = popenout.read()
  popenout.close()
  identity = float(out.split()[-1])
  if identity > min: return 1
  else: return 0
  

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  added = 0
  seqhash = defaultdict(str)
  ext = args['ext']
  for file in os.listdir('.'):
    if added == 1500: break
    if not file.endswith(ext): continue
    if not aln_is_conserved(file): continue
    fo = open(file)
    for line in fo:
      line = line.rstrip()
      if line.startswith(">"):
        id = line[1:]
        if id.count(" ") > 0: id = id[:id.index(" ")]
      else:
        seqhash[id] += line
    fo.close()
    added += 1
  for id, seq in seqhash.iteritems():
    print ">" + id
    print seq

# =============================================================================
args = handle_arguments()
main( args )

