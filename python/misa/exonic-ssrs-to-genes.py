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
  stdout( "usage: " + sys.argv[0] + " -e <path> -g <path> ")
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -e        all ssrs in exons" )
  stdout( " -g        parsed gff for all drosophilas" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hg:e:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-g': args['gff'] = value
    if key == '-e': args['ssr'] = value
    
  if not args.has_key('gff'):
    stderr( "parsed gff argument missing." )
    show_help()
  elif not file_exists( args.get('gff') ):
    stderr( "parsed gff does not exist." )
    show_help()
  
  if not args.has_key('ssr'):
    stderr( "ssr file argument missing." )
    show_help()
  elif not file_exists( args.get('ssr') ):
    stderr( "ssr file does not exist." )
    show_help()
  
  return args

# =============================================================================
class Gene():
  def __init__(self, line):
    cols = line.rstrip().split("\t")
    self.species = cols.pop(0)
    self.id = cols.pop(0)
    self.chr = cols.pop(0)
    self.start = cols.pop(0)
    self.stop = cols.pop(0)
    self.strand = cols.pop(0)
    self.loc = self.species + "|" + self.chr

  def pos_in_gene(self, pos):
    if int(pos) >= int(self.start) and int(pos) <= int(self.stop): return 1
    else: return 0

  
# =============================================================================
def get_gene_features(file):
  hash = defaultdict(list)
  fo = open(file)
  for line in fo:
    g = Gene(line)
    hash[g.loc].append(g)
  fo.close()
  return hash

# =============================================================================
def get_ssrs(file):
    hash = defaultdict(list)
    fo = open(file)
    for line in fo: 
      if line.startswith("ID\t"): continue
      m = MisaSSR(line)
      hash[m.geneid].append(m)
    fo.close()
    return hash


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  locssrs = get_ssrs(args['ssr'])
  locgenes = get_gene_features(args['gff'])
  for loc, genes in locgenes.iteritems():
    for gene in genes:
      for ssr in locssrs[loc]:
        if gene.pos_in_gene(ssr.startpos): print gene.id

# =============================================================================
args = handle_arguments()
main( args )

