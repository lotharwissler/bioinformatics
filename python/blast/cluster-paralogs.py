#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself
import hashlib

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        parsed blast.out, the first two columns being paralog pair IDs" )
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
		stderr( "parsed blast file file missing." )
		show_help()
	if not file_exists( args.get('file') ):
		stderr( "parsed blast file file does not exist." )
		show_help()
		
	return args

# =============================================================================
def get_all_pairs(ifile):
  parhash = {}
  fo = open( args.get('file'))
  for line in fo:
    id1, id2 = line.strip().split("\t")
    if not parhash.has_key(id1): parhash[id1] = []
    parhash[id1].append(id2)
    if not parhash.has_key(id2): parhash[id2] = []
    parhash[id2].append(id1)
  return parhash

# =============================================================================
def get_edges(ifile):
  hash = {}
  fo = open(ifile)
  for line in fo:
    ids = line.strip().split("\t")[0:2]
    ids.sort()
    hash[string.join(ids, ",")] = 1
  fo.close()
  return hash

# =============================================================================
def edges_in_cluster(clgids, edgehash):
  edgecount = 0
  for i in range(len(clgids)):
    for j in range(len(clgids)):
      if j <= i: continue
      gids = [clgids[i], clgids[j]]
      gids.sort()
      key = string.join(gids, ",")
      if edgehash.has_key(key): edgecount += 1
  return edgecount
    
# =============================================================================
# =============================================================================
def main( args ):
  parhash = {}
  nclusters = 0
  parhash = get_all_pairs(args.get('file'))
  nodes = parhash.keys()
  edgehash = get_edges(args['file'])
  while len(nodes) > 0:
    nclusters += 1
    first = nodes.pop()
    members = [first]
    check = []
    check.extend(parhash[first])
    while len(check) > 0:
      c = check.pop()
      if c in members: continue
      members.append(c)
      if parhash.has_key(c):
        check.extend(parhash[c])
        check = list(set(check))
    cid = hashlib.md5(string.join(members, '')).hexdigest()
    out = [cid, str(len(members)), str(edges_in_cluster(members, edgehash))] + members
    print string.join(out, "\t")
    for m in members: 
      if m in nodes: nodes.remove(m)
  print >> sys.stderr, "clusters: %s" % nclusters
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
