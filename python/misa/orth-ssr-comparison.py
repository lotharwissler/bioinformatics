#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import hashlib
from low import *  # custom functions, written by myself
from misa import MisaSSR
from collections import defaultdict


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -a <path> -b <path> -o <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -a        misa output of species 1 (Dmel)" )
  stdout( " -b        misa output of species 2 (Dxxx)" )
  stdout( " -o        flybase dmel ortholog report file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "ha:b:o:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-a': args['misa1'] = value
    if key == '-b': args['misa2'] = value
    if key == '-o': args['orth'] = value
    
  if not args.has_key('misa1'):
    stderr( "misa1 file argument missing." )
    show_help()
  elif not file_exists( args.get('misa1') ):
    stderr( "misa1 file does not exist." )
    show_help()

  if not args.has_key('misa2'):
    stderr( "misa2 file argument missing." )
    show_help()
  elif not file_exists( args.get('misa2') ):
    stderr( "misa2 file does not exist." )
    show_help()
 
  if not args.has_key('orth'):
    stderr( "orth file argument missing." )
    show_help()
  elif not file_exists( args.get('orth') ):
    stderr( "orth file does not exist." )
    show_help()
  
  return args


def get_orthologs(file, spec2):
  spec2 = spec2.lower()
  orthologs = {}
  fo = open(file)
  for line in fo:
    if line.startswith("#"): continue
    if len(line.rstrip()) == 0: continue
    columns = line.rstrip().split("\t")
    #print columns
    id1, id2, descr = columns[0], columns[5], columns[6]
    orthspecies = descr[:descr.index("\\")].lower()
    if orthspecies != spec2: continue
    orthologs[id1] = id2
  fo.close()
  return orthologs


def get_ssrs(file):
  hash = defaultdict(list)
  fo = open(file)
  for line in fo:
    if line.startswith("ID\t"): continue
    m = MisaSSR(line)
    hash[m.geneid].append(m)
  fo.close()
  return hash


def hash(s):
  return hashlib.sha224(s).hexdigest()

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  spec1 = args['misa1'][:args['misa1'].index("-")]
  spec2 = args['misa2'][:args['misa2'].index("-")]
  orthohash = get_orthologs(args['orth'], spec2)
  ssrs1 = get_ssrs(args['misa1'])
  ssrs2 = get_ssrs(args['misa2'])

  perfect, poly, shift, loss = 0, 0, 0, 0
  total = 0
  for gid1, ssrs in ssrs1.iteritems():
    if not orthohash.has_key(gid1): continue
    total += len(ssrs)
    gid2 = orthohash[gid1]
    if not ssrs2.has_key(gid2): 
      loss += len(ssrs)
      continue

    ossrs = ssrs2[gid2]

    # stage 1: perfect matches
    caught = {}
    for m1 in ssrs: 
      for m2 in ossrs:
        if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
        if m1.is_perfect_match_to(m2):
#          print "\nperfect match"
#          print m1.to_s()
#          print m2.to_s()
          perfect += 1
          caught[hash(m1.to_s())] = 1
          caught[hash(m2.to_s())] = 1

    # stage 2: polymorphic matches (same motif, but different number of repeats)
    for m1 in ssrs: 
      for m2 in ossrs:
        if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
        if m1.is_polymorphic_to(m2):
#          print "\npolymorphic match"
#          print m1.to_s()
#          print m2.to_s()
          poly += 1
          caught[hash(m1.to_s())] = 1
          caught[hash(m2.to_s())] = 1

    # stage 3: shifted matches (motif is shifted [permuated])
    for m1 in ssrs: 
      for m2 in ossrs:
        if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
        if m1.is_shifted_to(m2):
#          print "\nshifted match"
#          print m1.to_s()
#          print m2.to_s()
          shift += 1
          caught[hash(m1.to_s())] = 1
          caught[hash(m2.to_s())] = 1

    mapped = len(caught) / 2
#    print "\nUncaught between genes", gid1, "and", gid2
#    for m1 in ssrs:
#      if caught.has_key(hash(m1.to_s())): continue
#      print spec1 + "\t" + m1.to_s()
#
#    for m2 in ossrs:
#      if caught.has_key(hash(m2.to_s())): continue
#      print spec2 + "\t" + m2.to_s()
    loss += len(ssrs) - mapped 

  print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %( spec1, spec2, perfect, poly, shift, loss, total )




# =============================================================================
args = handle_arguments()
main( args )

