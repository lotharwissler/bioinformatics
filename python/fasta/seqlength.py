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
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file to import" )
  stdout( " -g        map file, tab delimited, regex to name (one per line) to group sequences into distinct bins" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:g:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-g': args['group'] = value
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
    
  return args


# =============================================================================
def read_groups( file ):
  groups = {}
  fo = open( file )
  for line in fo:
    line = line.rstrip()
    regex, name = line.split("\t")
    groups[name] = re.compile(regex)
  fo.close()
  return groups

# =============================================================================
def read_sequences( file, groups ):
  def add_entry( hash, groups, id, seq ):
    group = "*all*"
    for name, regex in groups.iteritems():
      if re.search(regex, id):
        group = name
        break
    if hash[group].has_key(id): sys.stderr.write("WARNING: overwriting entry with the same ID (%s) in group %s...\n" %(id, group))
    hash[group][id] = seq
    return hash


  hash = {}
  for name, regex in groups.iteritems(): hash[name] = {}
  if hash.has_key('*all*'): sys.stderr.write("WARNING: you used \"*all*\" as a group name. This name refers to all non-group-matching entries as well!\n")
  hash['*all*'] = {}

  id, seq = "", ""
  fo = open( file )
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"):
      if id != "": add_entry( hash, groups, id, seq )
      id = line[1:]
      seq = ""
    else:
      seq += line
  if id != "": add_entry( hash, groups, id, seq )
  fo.close()
  return hash

# =============================================================================
def eval_seq_lengths(hash):
  for group, seqhash in hash.iteritems():
    for id, seq in seqhash.iteritems():
      print string.join([group, id, str(len(seq))], "\t")

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  groups = {}
  if args.has_key('group'): groups = read_groups( args.get('group') )
  seqhash = read_sequences( args.get('file'), groups )
  eval_seq_lengths(seqhash)

# =============================================================================
args = handle_arguments()
main( args )

