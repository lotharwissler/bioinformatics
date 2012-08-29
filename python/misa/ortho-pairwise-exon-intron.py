#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import hashlib
from low import *  # custom functions, written by myself
from misa import MisaSSR
import newick
from collections import defaultdict
import pickle


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -a <path> -b <path> -o <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        combined misa output" )
  stdout( " -o        pairwise ortholog intra/intergenic regions file" )
  stdout( " -t        newick treew with branch lengths" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:o:t:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['misa'] = value
    if key == '-o': args['orth'] = value
    if key == '-t': args['tree'] = value
    
  if not args.has_key('misa'):
    stderr( "misa file argument missing." )
    show_help()
  elif not file_exists( args.get('misa') ):
    stderr( "misa file does not exist." )
    show_help()

  if not args.has_key('orth'):
    stderr( "orth file argument missing." )
    show_help()
  elif not file_exists( args.get('orth') ):
    stderr( "orth file does not exist." )
    show_help()
  
  return args


def get_distances(file):
  tree = open(file).readline().strip()
  ancestral_nodes = []
  leaves = {}
  while 1:
    # END OF TREE: semicolon
    if tree.startswith(";"): break

    # START INNER NODE
    if tree.startswith("("):
      tree = tree[1:]
      n = newick.Node()
      if len(ancestral_nodes) > 0: n.parent = ancestral_nodes[-1]
      ancestral_nodes.append(n)
      continue

    # END INNER NODE
    if tree.startswith(")"):
      tree = tree[1:]
      if re.match(":(\d+)", tree):
        distance = re.match(":(\d+)", tree).group(1)
        ancestral_nodes[-1].distance_to_parent = distance
        while re.match("[:\d]+", tree): tree = tree[1:]
      ancestral_nodes.pop(-1)
      continue

    # OUTER NODE SINGLE
    if re.match(",([A-Za-z]+):(\d+)\)", tree):
      els = re.match(",([A-Za-z]+):(\d+)", tree).groups()
      n1 = newick.Node()
      n1.parent = ancestral_nodes[-1]
      n1.distance_to_parent = els[1]
      leaves[els[0]] = n1
      while not tree.startswith(")"): tree = tree[1:]
      continue

    # OUTER NODE DOUBLE
    if re.match("([A-Za-z]+):(\d+),([A-Za-z]+):(\d+)", tree):
      els = re.match("([A-Za-z]+):(\d+),([A-Za-z]+):(\d+)", tree).groups()
      n1 = newick.Node()
      n1.parent = ancestral_nodes[-1]
      n1.distance_to_parent = els[1]
      n1.distance_to_parent = els[1]
      n2 = newick.Node()
      n2.parent = ancestral_nodes[-1]
      n2.distance_to_parent = els[3]
      leaves[els[0]] = n1
      leaves[els[2]] = n2
      while not tree.startswith(")"): tree = tree[1:]
      continue

    # INTERNAL INNER NODE
    if tree.startswith(",("):
      tree = tree[2:]
      n = newick.Node()
      if len(ancestral_nodes) > 0: n.parent = ancestral_nodes[-1]
      ancestral_nodes.append(n)
      continue
    if tree.startswith(","):
      tree = tree[1:]
      continue

  distances = {}
  for species1, leafnode1 in leaves.iteritems():
    for species2, leafnode2 in leaves.iteritems():
      distances[species1 + "," + species2] = str(leafnode1.summed_distance_to(leafnode2))
  return distances



class LocationPair():
  def __init__(self, line):
    columns = line.rstrip().split("\t")
    self.species = columns[0:2]
    self.type = columns[2]
    self.locations = [{'chr': columns[3], 'start': int(columns[4]), 'stop': int(columns[5])}, {'chr': columns[6], 'start': int(columns[7]), 'stop': int(columns[8])}]


def get_orthologs(file):
  orthologs = []
  fo = open(file)
  for line in fo:
    if line.startswith("#"): continue
    if len(line.rstrip()) == 0: continue
    orthologs.append(LocationPair(line))
  fo.close()
  return orthologs


def get_ssrs(file):
  hash = {}
  fo = open(file)
  for line in fo:
    if line.startswith("ID\t"): continue
    m = MisaSSR(line)
    hash[m.geneid + "|" + str(m.startpos)] = m
  fo.close()
  return hash


def hash(s):
  return hashlib.sha224(s).hexdigest()

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  ssrs = get_ssrs(args['misa'])
  orthologLocationPairs = get_orthologs(args['orth'])
  distances = get_distances(args['tree'])

  perfect, poly, shift, loss = defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int)), defaultdict(lambda: defaultdict(int))
  for locpair in orthologLocationPairs:
    qspecies, ospecies = locpair.species[0], locpair.species[1]
    qchr, qstart, qstop = locpair.locations[0]['chr'], locpair.locations[0]['start'], locpair.locations[0]['stop']
    ochr, ostart, ostop = locpair.locations[1]['chr'], locpair.locations[1]['start'], locpair.locations[1]['stop']
    qssrs, ossrs = [], []
    for s in range(qstart, qstop):
      key = qspecies + "|" + qchr + "|" + str(s)
      if ssrs.has_key(key): qssrs.append(ssrs[key])
    for s in range(ostart, ostop):
      key = ospecies + "|" + ochr + "|" + str(s)
      if ssrs.has_key(key): ossrs.append(ssrs[key])
    key = [qspecies, ospecies]
    key.sort()
    key = string.join(key, ",")

    # no SSRs in these both locations
    if len(qssrs) == 0 and len(ossrs) == 0: continue
    # no SSRs in either one of the two locations
    if len(qssrs) == 0:
      for ssr in ossrs: loss[key][ssr.feature] += 1
      continue
    if len(ossrs) == 0:
      for ssr in qssrs: loss[key][ssr.feature] += 1
      continue

    caught = {}
    # stage 1: perfect matches
    for m1 in qssrs: 
      for m2 in ossrs:
        if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
        if m1.is_perfect_match_to(m2):
          if m1.feature == m2.feature: 
            perfect[key][m1.feature] += 1
            if m1.feature != "intergenic": 
              print m1.to_s()
              print m2.to_s()
            caught[hash(m1.to_s())] = 1
            caught[hash(m2.to_s())] = 1

    # stage 2: polymorphic matches (same motif, but different number of repeats)
    for m1 in qssrs: 
      for m2 in ossrs:
        if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
        if m1.is_polymorphic_to(m2):
          if m1.feature == m2.feature: 
            poly[key][m1.feature] += 1
            if m1.feature != "intergenic": 
              print m1.to_s()
              print m2.to_s()
            caught[hash(m1.to_s())] = 1
            caught[hash(m2.to_s())] = 1

    # stage 3: shifted matches (motif is shifted [permuated])
    for m1 in qssrs: 
      for m2 in ossrs:
        if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
        if m1.is_shifted_to(m2):
          if m1.feature == m2.feature: 
            shift[key][m1.feature] += 1
            if m1.feature != "intergenic": 
              print m1.to_s()
              print m2.to_s()
            caught[hash(m1.to_s())] = 1
            caught[hash(m2.to_s())] = 1

    for m1 in qssrs:
      if caught.has_key(hash(m1.to_s())): continue
      loss[key][m1.feature] += 1
    for m2 in ossrs:
      if caught.has_key(hash(m2.to_s())): continue
      loss[key][m2.feature] += 1 
  
  
  keys = list(set(perfect.keys() + poly.keys() + shift.keys() + loss.keys()))
  keys.sort()
  for key in keys:
    speciespair = key
    time = str(distances[speciespair])
    perfectcounts, polycounts, shiftedcounts, losscounts = [], [], [], []
    for feature in ["exon", "intron", "intergenic"]:
      count = "0"
      if loss[key].has_key(feature): count = str(loss[key][feature])
      losscounts.append(count)
      
      count = "0"
      if perfect[key].has_key(feature): count = str(perfect[key][feature])
      perfectcounts.append(count)

      count = "0"
      if poly[key].has_key(feature): count = str(poly[key][feature])
      polycounts.append(count)

      count = "0"
      if shift[key].has_key(feature): count = str(shift[key][feature])
      shiftedcounts.append(count)

    sys.stderr.write(string.join([key, speciespair, time], "\t"))
    sys.stderr.write("\t" + string.join(perfectcounts, "\t"))
    sys.stderr.write("\t" + string.join(polycounts, "\t"))
    sys.stderr.write("\t" + string.join(shiftedcounts, "\t"))
    sys.stderr.write("\t" + string.join(losscounts, "\t") + "\n")


# =============================================================================
args = handle_arguments()
main( args )

