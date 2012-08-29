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
  stdout( "usage: " + sys.argv[0] + " -f <path> -g <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        GO flat file to import [tab delimited]" )
  stdout( " -g        gene_id to group_id table [tab delimited]" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:g:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-g': args['group'] = value
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
    
  if not args.has_key('group'):
    stderr( "group file argument missing." )
    show_help()
  elif not file_exists( args.get('group') ):
    stderr( "group file does not exist." )
    show_help()
 
  return args


# =============================================================================
def get_gene2groups(file):
  hash = {}
  groups = {}
  fo = open(file)
  for line in fo:
    if line.startswith("#"): continue
    if not len(line.split("\t")) == 2: continue
    geneid, group = line.rstrip().split("\t")
    if not hash.has_key(geneid): hash[geneid] = []
    hash[geneid].append(group)
    groups[group] = 1
  fo.close()
  return hash, groups.keys()

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  gene2groups, groups = get_gene2groups(args['group'])
  groups.sort()

  hash = {}
  fo = open( args.get('file') )
  for line in fo:
    line = line.strip()
    geneid, goterm = line.split("\t")
    if geneid.count(" ") > 0: geneid = geneid[:geneid.index(" ")]
    if not gene2groups.has_key(geneid): continue
    if not hash.has_key(goterm): hash[goterm] = defaultdict(int)
    for g in gene2groups[geneid]: hash[goterm][g] += 1
  fo.close()
  print string.join(["GO.term"] + groups, "\t")
  for goterm, counthash in hash.iteritems():
    #print goterm
    #print hash[goterm]
    #for g in groups:
    #  print groups, g, counthash[g]
    counts = [counthash[g] for g in groups]
    if sum(counts) < 5: continue
    counts = [str(c) for c in counts]
    print string.join([goterm] + counts, "\t")

# =============================================================================
args = handle_arguments()
main( args )

