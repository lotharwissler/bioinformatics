#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from collections import defaultdict
from misa import MisaSSRspecies


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        all.misa out file" )
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
    stderr( "fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "fasta file does not exist." )
    show_help()
  
  return args

  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  specieshash = {}
  fo = open(args['file'])
  for line in fo:
    m = MisaSSRspecies(line)
    if not specieshash.has_key(m.species): specieshash[m.species] = {'gc':defaultdict(int), 'p1':0, 'p2':0, 'p3':0, 'p4':0, 'p5':0, 'p6':0, 'p1l':0, 'p2l':0, 'p3l':0, 'p4l':0, 'p5l':0, 'p6l':0, 'ssrl':0}
    # gc
    for char in ['A', 'T', 'G', 'C']:
      specieshash[m.species]['gc'][char] += m.motif.count(char) * m.repeats
    # count repeats and coverage
    specieshash[m.species][m.type] += 1
    specieshash[m.species][m.type + 'l'] += m.length
    specieshash[m.species]['ssrl'] += m.length

  speciesarray = specieshash.keys()
  speciesarray.sort()
  print "#species\tssr.gc\tp1\tp2\tp3\tp4\tp5\tp6\tp1l\tp2l\tp3l\tp4l\tp5l\tp6l\tp1a\tp2a\tp3a\tp4a\tp5a\tp6a\tssrl"
  for species in speciesarray:
    total = sum(specieshash[species]['gc'].values())
    gc = 1.0 * (specieshash[species]['gc']['G'] + specieshash[species]['gc']['C']) / total
    repeats = [specieshash[species]['p1'], specieshash[species]['p2'], specieshash[species]['p3'], specieshash[species]['p4'], specieshash[species]['p5'], specieshash[species]['p6']]
    repeats = [str(r) for r in repeats]
    coverage = [specieshash[species]['p1l'], specieshash[species]['p2l'], specieshash[species]['p3l'], specieshash[species]['p4l'], specieshash[species]['p5l'], specieshash[species]['p6l']]
    coverage = [str(c) for c in coverage]
    avglength = []
    for i in range(len(repeats)): avglength.append(str(float(coverage[i]) / float(repeats[i])))
    print species + "\t" + str(gc) + "\t" + string.join(repeats, "\t") + "\t" + string.join(coverage, "\t") + "\t" + string.join(avglength, "\t") + "\t" + str(specieshash[species]['ssrl'])

# =============================================================================
args = handle_arguments()
main( args )

