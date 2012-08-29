#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
import gff3
import fasta
from Bio.Seq import Seq

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -i        surrounding genes (ID1,ID2)" )
  stdout( " -g        gff file" )
  stdout( " -f        fasta file" )
  stdout( " -r        return reverse complement" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """

  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hi:g:f:r" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'rv': False}
  for key, value in keys:
    if key == '-i': args['genes'] = value.split(',')
    if key == '-g': args['gffile'] = value
    if key == '-f': args['fastafile'] = value
    if key == '-r': args['rv'] = True
    
  for key in ['genes', 'gffile', 'fastafile']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key):
      print >> sys.stderr, "missing argument", key
      show_help()
  return args

# =============================================================================
def statusbar(current, total, message="", width=40):
  progress = 1.0*current/total
  if message != "": message = "[" + message + "]"
  progressbar = "=" * int(progress*width)
  while len(progressbar) < width: progressbar += " " 
  sys.stderr.write("\r   0% " + progressbar + " 100% " + message)
  if progress == 1.0: sys.stderr.write("\n")
  

# =============================================================================
def get_coordinates(gfhash, genes):
  positions = []
  for scaffold, gfs in gfhash.iteritems():
    for gf in gfs:
      if not gf.get_attributes().has_key('ID'): continue
      if gf.get_attributes()['ID'] in genes: positions += [gf.start, gf.stop]
      if len(positions) == 4: break
    if len(positions) == 4: break
  if min(positions[0:2]) < min(positions[2:4]):
    return scaffold, max(positions[0:2]), min(positions[2:4])
  else:
    return scaffold, max(positions[2:4]), min(positions[0:2])
        
    
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  gfhash = gff3.get_gff_hash(args['gffile'])
  sys.stderr.write("gff loaded ")
  gid, startpos, endpos = get_coordinates(gfhash, args['genes'])
  sys.stderr.write("| coordinates identified ")
  if not args['rv']: print ">%s_%s:%s" %(gid, startpos, endpos)
  else: print ">%s_%s:%s" %(gid, endpos, startpos)
  
  seqhash = fasta.get_sequence_hash(args['fastafile'])
  sys.stderr.write("| fasta loaded ")
  seq = seqhash[gid][startpos-1:endpos]
  if args['rv']: seq = Seq(seq).reverse_complement().tostring()
  sys.stderr.write("| subsequence extracted ")
  print seq
  sys.stderr.write("\n")

# =============================================================================
args = handle_arguments()
main( args )

