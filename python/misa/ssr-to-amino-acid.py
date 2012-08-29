#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself
from misa import MisaSSR

AA = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print >> sys.stderr, "usage: " + sys.argv[0] + " -d <gff-folder>"
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -p        protein fasta file" )
  stdout( " -m        misa file incl. protein in last column" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hm:p:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-m': args['misa'] = value
    if key == '-p': args['protein'] = value
    
  if not args.has_key('misa'):
    print >> sys.stderr, "misa file argument missing."
    show_help()
  elif not file_exists( args.get('misa') ):
    print >> sys.stderr, "misa file does not exist."
    show_help()

  if not args.has_key('protein'):
    print >> sys.stderr, "protein file argument missing."
    show_help()
  elif not file_exists( args.get('protein') ):
    print >> sys.stderr, "protein file does not exist."
    show_help()

  return args


# =============================================================================
def get_ssrs(file):
    hash = defaultdict(list)
    fo = open(file)
    for line in fo: 
      if line.startswith("ID\t"): continue
      m = MisaSSR(line)
      hash[m.geneid].append(m)
    fo.close()
    print >> sys.stderr, "read %s microsatellites" % len(hash)
    return hash

# =============================================================================
def get_protein(file):
  seqhash = defaultdict(str)
  id = ""
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"):
      id = line[1:]
      if id.count(" ") > 0: id = id[:id.index(" ")]
    else:
      seqhash[id] += line
  fo.close()
  print >> sys.stderr, "read %s protein sequencess" % len(seqhash)
  return seqhash

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  SSRs = get_ssrs(args['misa']) 
  seqhash = get_protein(args['protein']) 
  for sid, SSRs in SSRs.iteritems():
    for SSR in SSRs:
      prot = seqhash[SSR.feature]
      pstart, pend = SSR.startpos / 3, SSR.endpos / 3
      seq = prot[pstart:pend+1]
      indic = "-"
      for aa in AA:
        if seq.count(aa*4) > 0:
          indic = aa
          break
      if indic == "-":
        for aa1 in AA:
          for aa2 in AA:
            if seq.count((aa1+aa2)*3) > 0:
              indic = aa1+aa2
              break
      print SSR.to_s() + "\t" + indic + "\t" + seq


# =============================================================================
args = handle_arguments()
main( args )
