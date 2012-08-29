#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -n <file> -p <file> [-i <string>]" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -n        nucleotide fasta file (trancript, CDS, ...)" )
  stdout( " -p        protein fasta file" )
  stdout( " -i        string to prefix all unique IDs with" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hp:n:i:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  args['idpref'] = ""
  for key, value in keys:
    if key == '-p': args['pfasta'] = value
    if key == '-n': args['nfasta'] = value
    if key == '-i': args['idpref'] = value

    
  if not args.has_key('pfasta'):
    stderr( "protein fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('pfasta') ):
    stderr( "protein fasta file does not exist." )
    show_help()

  if not args.has_key('nfasta'):
    stderr( "nucleotide fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('nfasta') ):
    stderr( "nucleotide fasta file does not exist." )
    show_help()

  return args


def get_all_ids_from_fasta(file):
  ids = []
  fo = open(file)
  for line in fo:
    line = line.strip()
    if not line.startswith(">"): continue
    ids.append(line[1:])
  return ids


def get_unique_id_fragments(idArray):
  hash = defaultdict(int)
  part2id = {}
  for id in idArray:
    parts = id.split("|")
    for part in parts: 
      hash[part] += 1
      part2id[part] = id
  
  uniqHash = {}
  for part, count in hash.iteritems():
    if count == 1: uniqHash[part] = part2id[part]
  revUniqHash = defaultdict(list)
  
  for part, id in uniqHash.iteritems():
    revUniqHash[id].append(part)
  return uniqHash, revUniqHash


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  allPids = get_all_ids_from_fasta(args.get('pfasta'))
  allNids = get_all_ids_from_fasta(args.get('nfasta'))
  upPids, uHashPids = get_unique_id_fragments(allPids)
  upNids, uHashNids = get_unique_id_fragments(allNids)
  
  matchHash = {}
  matchPart = {}
  prefix = args.get('idpref')
  for pid, uparts in uHashPids.iteritems():
    debug = 0
    if pid == "jgi|Araly1|878105|Al_scaffold_0002_2699": debug = 1
    for upart in uparts:
      if debug: print >> sys.stderr, upart, pid
      if upNids.has_key(upart):
        nid = upNids[upart]
        if debug: print >> sys.stderr, "match to", nid
        if matchHash.has_key(pid):
          if debug: print >> sys.stderr, "ERROR: matching Hash already contains an association with PID", pid, "=>", matchHash[pid]
          continue
        if matchHash.has_key(nid):
          if debug: print >> sys.stderr, "ERROR: matching Hash already contains an association with NID", nid, "=>", matchHash[nid]
          continue
        matchHash[pid] = nid
        matchHash[nid] = pid
        matchPart[pid + "$$$" + nid] = upart
        print string.join([prefix+upart, pid, nid], "\t")
        break
    if not matchHash.has_key(pid):
      print >> sys.stderr, "no match for PID", pid

  if len(matchPart) == len(allPids) and len(matchPart) == len(allNids):
    print >> sys.stderr, "everything is well: we have found an association match between all protein and nucleotide IDs"
  else:
    print >> sys.stderr, "ERROR: unequal number of matches (%s) and protein/nucleotide IDs (%s/%s)" %( len(matchPart), len(allPids), len(allNids) )

  
# =============================================================================
args = handle_arguments()
main( args )

