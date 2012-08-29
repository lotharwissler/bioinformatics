#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
import fasta
import numpy
import copy

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file (scaffolds/chromosomes)" )
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
    if key == '-f': args['fastafile'] = value
    
  for key in ['fastafile']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key):
      print >> sys.stderr, "missing argument", key
      show_help()
  return args

# =============================================================================
def statusbar(current, total, message="", width=40):
  progress = 1.0*current/total
  if message != "": message = "[" + message + "]"
  progressbar = "=" * int(progress*width)
  while len(progressbar) < width: progressbar += " " 
  sys.stderr.write("\r   0% " + progressbar + " 100% " + message)
  if progress == 1.0: sys.stderr.write("\n")
  
# =============================================================================
def get_n50_from_lengthhash(sid2length):
  sorted_lengths = sid2length.values()
  sorted_lengths.sort()
  total = sum(sid2length.values())
  i, runningsum, threshold = -1, 0, 0.5*total
  while runningsum < threshold:
    i += 1
    runningsum += sorted_lengths[i]
  return sorted_lengths[i-1]

# =============================================================================
def get_gaps_and_Ns(fastafile):
  gaps, n = 0, 0
  fo = open(fastafile)
  for line in fo:
    if line.startswith(">"): continue
    gaps += line.count("-")
    n += line.count("N")
  fo.close()
  return gaps, n
  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  sid2length = fasta.get_length_hash(args['fastafile'])
  count = len(sid2length)
  lmin, lmax = min(sid2length.values()), max(sid2length.values())
  mean, median = numpy.mean(sid2length.values()), numpy.median(sid2length.values())
  n50 = get_n50_from_lengthhash(sid2length)
  total = sum(sid2length.values())
  gaps, unresolved = get_gaps_and_Ns(args['fastafile'])
  print string.join([str(e) for e in [args['fastafile'], total, count, lmin, lmax, mean, median, n50, gaps, unresolved]], "\t")

# =============================================================================
args = handle_arguments()
main( args )

