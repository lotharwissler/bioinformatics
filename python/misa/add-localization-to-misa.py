#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself
from misa import MisaSSRspecies
import pickle

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print >> sys.stderr, "usage: " + sys.argv[0] + " -d <gff-folder>"
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -g        gff3 file" )
  stdout( " -f        misa file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hg:f:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-g': args['gff'] = value
    if key == '-f': args['misa'] = value
    
  if not args.has_key('gff'):
    print >> sys.stderr, "gff file argument missing."
    show_help()
  elif not file_exists( args.get('gff') ):
    print >> sys.stderr, "gff file does not exist."
    show_help()

  if not args.has_key('misa'):
    print >> sys.stderr, "misa file argument missing."
    show_help()
  elif not file_exists( args.get('misa') ):
    print >> sys.stderr, "misa file does not exist."
    show_help()

  return args

# =============================================================================
def get_ssrs(file):
    hash = defaultdict(list)
    fo = open(file)
    for line in fo: 
      if line.startswith("ID\t"): continue
      m = MisaSSRspecies(line)
      hash[m.species + '|' + m.geneid].append(m)
    fo.close()
    return hash

# =============================================================================
def get_features(filename):
  type2abbrv = { 'exon':'E', 'intron':'I', 'five_prime_UTR':'5', 'three_prime_UTR':'3' }
  features = {}
  fo = open(filename)
  for line in fo: 
    if line.startswith("#") or len(line.rstrip()) == 0: continue
    columns = line.rstrip().split("\t")
    if len(columns) != 9: continue
    type = columns[2]
    if type != "sequence_assembly" and type != "exon" and type != "intron" and type != "five_prime_UTR" and type != "three_prime_UTR": continue
    chr, start, stop, strand, descr = columns[0], columns[3], columns[4], columns[6], columns[8]
    if type == "sequence_assembly": features[chr] = ["i"] * int(stop)
    else:
      for i in range(int(start)-1, int(stop)-1): features[chr][i] = type2abbrv[type]
  fo.close()
  print >> sys.stderr, "features of all %s scaffolds loaded." % len(features)
  return features

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  Features = get_features(args['gff']) 
  fo = open(args['misa'])
  for line in fo:
    if line.startswith("ID\t"): continue
    line = line.rstrip()
    columns = line.split("\t")
    key = columns[1]
    start, stop = int(columns[6])-1, int(columns[7])-1
    fstart, fstop = Features[key][start], Features[key][stop]
    if fstart != fstop: print >> sys.stderr, "SSR spans two different features: %s %s/%s" %( key, fstart, fstop )
    print line + "\t" + fstart
  fo.close()
  


# =============================================================================
args = handle_arguments()
main( args )
