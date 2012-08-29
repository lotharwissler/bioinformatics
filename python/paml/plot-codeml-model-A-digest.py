#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import math        # match functions
from low import *  # custom functions, written by myself
from collections import defaultdict

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        digested file with model A positional data of PS" )
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
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
  
  return args


# =============================================================================
def get_genes_with_features( file ):
  hash = defaultdict(list)
  fo = open(file)
  for line in fo:
    filename, length, feature = line.rstrip().split("\t")[0:3]
    cluster = filename[:filename.index(".")]
    branch = filename[filename.index("tree.")+5:filename.rindex(".")]
    pos, aa, prob = feature.split()[0:3]
    if prob.endswith("*"): prob = prob[:prob.index("*")]
    prob = float(prob)
    pos = int(pos)
    length = int(length)
    hash[cluster + " " + branch].append([length, pos, aa, prob])
  fo.close()
  return hash


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  hash = get_genes_with_features(args['file'])
  for key, featurearray in hash.iteritems():
    cluster, branch = key.split()
    length = int(featurearray[0][0])
    import matplotlib.pyplot as P
    x = [e+1 for e in range(length+1)]
    y1 = [0] * (length+1)
    y2 = [0] * (length+1)
    for feature in featurearray:
      length, pos, aa, prob = feature[0:4]
      if prob > 0.95: y1[pos] = prob
      else: y2[pos] = prob
    
    P.bar(x, y1, color='#000000', edgecolor='#000000')
    P.bar(x, y2, color='#bbbbbb', edgecolor='#bbbbbb')
    P.ylim(ymin=0, ymax=1)
    P.xlim(xmin=0, xmax=length)
    P.xlabel("position in the ungapped alignment [aa]")
    P.ylabel(r'$P (\omega > 1)$')
    P.title(cluster + " (branch " + branch + ")")

    P.axhline(y=.95, xmin=0, xmax=length, linestyle=":", color="k")
    P.savefig(cluster + "." + branch + ".png", format="png")
    P.close()


# =============================================================================
args = handle_arguments()
main( args )

