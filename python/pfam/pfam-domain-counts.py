#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from collections import defaultdict
import glob
import pfam
import stats


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -h <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        hmmout dir (*.hmmout)" )
  stdout( " -p        protein ids dir (*.gids; optional, if given only domains in the specified proteins are considered)" )
  stdout( " -s        list of species (csv)" )
  stdout( " -i        pfam dids to ignore" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """

  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "h:p:s:i:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-h': args['hmmoutdir'] = value
    if key == '-i': args['ignorefile'] = value
    if key == '-p': args['pidsdir'] = value
    if key == '-s': args['species'] = value.split(',')
    
  for key in ['hmmoutdir', 'species']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key):
      print >> sys.stderr, "missing argument", key
      show_help()
  return args

  
# =============================================================================
def get_pids(ifile):
  ohash = {}
  fo = open(ifile)
  for line in fo:
    gid = line.rstrip()
    ohash[gid] = 1
  fo.close()
  return ohash



# =============================================================================
def did2count(hmmout, pids=False):
  pid2pfamdomains = pfam.read_hmmout(hmmout)
  if not pids: keys = pid2pfamdomains.keys()
  else: keys = pids
  did2count = {}
  for pid in keys:
    if not pid2pfamdomains.has_key(pid): continue
    for did in list(set([d.get_attr("hmm_name") for d in pid2pfamdomains[pid]])):
      if not did2count.has_key(did): did2count[did] = 0
      did2count[did] += 1
  return did2count
  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  species2counts = {}
  for species in args['species']:
    hmmoutfiles = glob.glob(args['hmmoutdir'] + "/" + species + "*.hmmout")
    if not len(hmmoutfiles) == 1: sys.exit("ERROR: no single hmmout file found for species %s" % species)
    if args.has_key('pidsfile'):  
      pidsfiles = glob.glob(args['pidsdir'] + "/" + species + "*.gids")
      if not len(pidsfiles) == 1: sys.exit("ERROR: no single pids file found for species %s" % species)
      pids = get_pids(pidsfiles[0]).keys()
    else: pids = False
    species2counts[species] = did2count(hmmoutfiles[0], pids)
  
  ignore = {}
  if args.has_key('ignorefile'): ignore = get_pids(args['ignorefile'])
  hash = {}
  for s, counts in species2counts.iteritems():
    for did in counts.keys(): 
      if not ignore.has_key(did): hash[did] = 1
    
  fwt = open("pfam-table", "w")
  fwt.write(string.join(["DID"] + args['species'], "\t") + "\n")
  for did in hash.keys():
    out = did
    for species in args['species']:
      count = 0
      if species2counts[species].has_key(did): count = species2counts[species][did]
      out += "\t" + str(count)
    fwt.write(out + "\n")
  fwt.close()
  
  fwm = open("pfam-matrix.csv", "w")
  fwm.write("," + string.join(args['species'], ",") + "\n")
  for i in range(len(args['species'])):
    s1 = args['species'][i]
    fwm.write(s1)
    for j in range(len(args['species'])):
      if i == j: 
        fwm.write(",1")
        continue
      s2 = args['species'][j]
      v1, v2 = [], []
      for did in hash.keys():
        v1.append(species2counts[s1].get(did, 0))
        v2.append(species2counts[s2].get(did, 0))
      cor, p = stats.correlate(v1, v2)
      fwm.write("," + str(cor))
      #print string.join([s1, s2, str(cor), str(p)], "\t")#
    fwm.write("\n")
  fwm.close()

# =============================================================================
args = handle_arguments()
main( args )

