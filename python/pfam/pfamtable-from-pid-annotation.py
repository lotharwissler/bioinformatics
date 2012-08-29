#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
import glob
import stats


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -d        pfam dir" )
  stdout( " -e        file extention (e.g. \"*.pfam\")" )
  stdout( " -m        output matrix instead of table" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """

  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hd:e:m" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'matrix': False}
  for key, value in keys:
    if key == '-d': args['pfamdir'] = value
    if key == '-e': args['extension'] = value
    if key == '-m': args['matrix'] = True
    
  for key in ['pfamdir', 'extension']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key):
      print >> sys.stderr, "missing argument", key
      show_help()
  return args

  
# =============================================================================
def get_domain_counts(ifile):
  chash = {}
  fo = open(ifile)
  for line in fo:
    pid, did = line.strip().split("\t")
    if not chash.has_key(did): chash[did] = 0
    chash[did] += 1
  fo.close()
  return chash


# =============================================================================
def gather_input(idir, ext):
  hash = {}
  for filename in glob.glob(idir + '/' + ext):
    species = os.path.split(filename)[1]
    species = species[:species.index(ext[1:])]
    hash[species] = filename
  return hash

# =============================================================================
def output_table(species2domain2count):
  species = species2domain2count.keys()
  species.sort()
  alldids = {}
  for spec, chash in species2domain2count.iteritems():
    for did, count in chash.iteritems(): alldids[did] = 1
    
  print string.join(["DID"] + species, "\t")
  for did in alldids.keys():
    out = did
    for spec in species: out += "\t" + str(species2domain2count[spec].get(did,0))
    print out
    

# =============================================================================
def score_pair(v1, v2, method=1):
  if method == 1:
    cor, p = stats.correlate(v1, v2)
    return cor
  elif method == 2:
    x1, x2 = [], []
    for i in range(len(v1)):
      if v1[i] != 0 or v2[i] != 0: x1.append(v1[i]), x2.append(v2[i])
    cor, p = stats.correlate(x1, x2)
    return cor
  elif method == 3:
    up = 0
    for i in range(len(v1)):
      if v1[i] != 0 and v2[i] != 0: up += 1
    return 1.0*up/len(v1)
  elif method == 4:
    up = 0
    for i in range(len(v1)):
      if v1[i] != 0 and v2[i] != 0: up += 1
    return up

# =============================================================================
def output_matrix(species2domain2count):
  species = species2domain2count.keys()
  species.sort()
  alldids = {}
  for spec, chash in species2domain2count.iteritems():
    for did, count in chash.iteritems(): alldids[did] = 1
    
  print string.join([""] + species, ",")
  for sp1 in species:
    out = sp1
    for sp2 in species:
      v1, v2 = [], []
      for did in alldids.keys():
        v1.append(species2domain2count[sp1].get(did,0))
        v2.append(species2domain2count[sp2].get(did,0))
      sim = score_pair(v1, v2, method=4)
      out += "," + str(sim)
    print out
      
     

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  species2file = gather_input(args['pfamdir'], args['extension'])
  species2domain2count = {}
  for species, ifile in species2file.iteritems(): species2domain2count[species] = get_domain_counts(ifile)
  alldids = {}
  for species, chash in species2domain2count.iteritems():
    for did, count in chash.iteritems(): alldids[did] = 1
  
  if args['matrix']: output_matrix(species2domain2count)
  else: output_table(species2domain2count)
  
  

# =============================================================================
args = handle_arguments()
main( args )

