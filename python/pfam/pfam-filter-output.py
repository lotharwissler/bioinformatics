#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import math        # match functions
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> " )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        pfam_full file to parse" )
  stdout( " -d        Pfam-A.hmm file [default: /global/databases/pfam/current/pfam_scan_db/Pfam-A.hmm" )
  stdout( " -c        cutoff to apply (GA|TC|NC) [default: GA]" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:d:c:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'hmmfile':'/global/databases/pfam/current/pfam_scan_db/Pfam-A.hmm', 'cut':'GA'}
  for key, value in keys:
    if key == '-f': args['annotfile'] = value
    if key == '-d': args['hmmfile'] = value
    if key == '-c': args['cut'] = value
    
  for key in ['annotfile', 'hmmfile']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key):
      print >> sys.stderr, "missing argument", key
      show_help()

  return args


# =============================================================================
def get_regex( args ):
  idhash = {}
  idhash['name'] = re.compile('^#=GF ID\s+(\S+)')
  idhash['acc'] = re.compile('^#=GF AC\s+(PF\S+)')
  idhash['descr'] = re.compile('^#=GF DE\s+(.*)$')
  idhash['comment'] = re.compile('^#=GF CC\s+(.*)$')
  idhash['pftype'] = re.compile('^#=GF TP\s+(\S+)')
  idhash['terminate'] = re.compile('^\\\\$')
  return idhash

# =============================================================================
class PfamEntry:
  def __init__(self):
    self.name = None
    self.acc = None
    self.descr = None
    self.ga = None
    self.tc = None
    self.nc = None
 

# =============================================================================
def load_hmmfile(hmmfile):
  name2pfam = {}
  fo = open( hmmfile )
  entry = PfamEntry()
  for line in fo:
    line = line.rstrip()
    if line.startswith('//'):
      name2pfam[entry.name] = entry
      entry = PfamEntry()
    elif line.startswith("NAME"): entry.name = line.split()[1]
    elif line.startswith("ACC"): entry.acc = line.split()[1]
    elif line.startswith("DESC"): entry.descr = line.split(' ', 1)[1]
    elif line.startswith("GA"): entry.ga = float(line.split()[1])
    elif line.startswith("NC"): entry.nc = float(line.split()[1])
    elif line.startswith("TC"): entry.tc = float(line.split()[1])
  fo.close()
  return name2pfam
  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  name2entry = load_hmmfile(args['hmmfile'])
  fo = open(args['annotfile'])
  for line in fo:
    line = line.rstrip()
    if line.startswith("#") or len(line) == 0: continue
    cols = line.split()
    name, score = cols[-9], float(cols[-4])
    if args['cut'] == 'GA' and score < name2entry[name].ga: continue
    elif args['cut'] == 'TC' and score < name2entry[name].tc: continue
    elif args['cut'] == 'NC' and score < name2entry[name].nc: continue
    print string.join(cols, "\t")
  fo.close()    

# =============================================================================
args = handle_arguments()
main( args )
