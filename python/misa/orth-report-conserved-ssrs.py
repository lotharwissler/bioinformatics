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


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -a <path> -b <path> -o <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        combined misa output" )
  stdout( " -o        clustered flybase dmel ortholog file" )
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
#      print "end of an"
      if re.match(":(\d+)", tree):
        distance = re.match(":(\d+)", tree).group(1)
        ancestral_nodes[-1].distance_to_parent = distance
#        print "   ... distance:", distance
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
#      print "single node found:", els[0], "distance:", els[1]
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
#      print "double node found:", els[0], els[2], "distances:", els[1], els[3]
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
      #ancestral_nodes[-1].parent = ancestral_nodes[-2]
      #ancestral_nodes.pop()
      tree = tree[1:]
      continue

  distances = {}
  for species1, leafnode1 in leaves.iteritems():
    for species2, leafnode2 in leaves.iteritems():
      distances[species1 + "," + species2] = str(leafnode1.summed_distance_to(leafnode2))
  return distances






def get_orthologs(file):
  orthologs = []
  fo = open(file)
  for line in fo:
    if line.startswith("#"): continue
    if len(line.rstrip()) == 0: continue
    columns = line.rstrip().split("\t")
    orthologs.append(columns)
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
  
  ssrs = get_ssrs(args['misa'])
  orthologclusters = get_orthologs(args['orth'])
  distances = get_distances(args['tree'])
  
  for cluster in orthologclusters:
    for i in range(len(cluster)):
      query = cluster[i]
      qid = query[:query.index("(")]
      qspecies = query[query.index("(")+1:-1]
      q_ssrs = ssrs[qid]
      for j in range(i+1,len(cluster)):
        ortholog = cluster[j]
        oid = ortholog[:ortholog.index("(")]
        ospecies = ortholog[ortholog.index("(")+1:-1]
        o_ssrs = ssrs[oid]
        key = [qspecies, ospecies]
        key.sort()
        key = string.join(key, ",")
        if len(q_ssrs) == 0 and len(o_ssrs) == 0: 
          break
        if len(q_ssrs) == 0:
          break
        if len(o_ssrs) == 0:
          break

        # stage 1: perfect matches
        caught = {}
        for m1 in q_ssrs: 
          for m2 in o_ssrs:
            if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
            if m1.is_perfect_match_to(m2):
#          print "\nperfect match"
#          print m1.to_s()
#          print m2.to_s()
              print "perfect\t" + m1.to_s()
              print "perfect\t" + m2.to_s()
              caught[hash(m1.to_s())] = 1
              caught[hash(m2.to_s())] = 1

        # stage 2: polymorphic matches (same motif, but different number of repeats)
        for m1 in q_ssrs: 
          for m2 in o_ssrs:
            if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
            if m1.is_polymorphic_to(m2):
#          print "\npolymorphic match"
#          print m1.to_s()
#          print m2.to_s()
              print "poly\t" + m1.to_s()
              print "poly\t" + m2.to_s()
              caught[hash(m1.to_s())] = 1
              caught[hash(m2.to_s())] = 1

        # stage 3: shifted matches (motif is shifted [permuated])
        for m1 in q_ssrs: 
          for m2 in o_ssrs:
            if caught.has_key(hash(m1.to_s())) or caught.has_key(hash(m2.to_s())): continue
            if m1.is_shifted_to(m2):
#          print "\nshifted match"
#          print m1.to_s()
#          print m2.to_s()
              print m1.to_s()
              print m2.to_s()
              caught[hash(m1.to_s())] = 1
              caught[hash(m2.to_s())] = 1



# =============================================================================
args = handle_arguments()
main( args )

