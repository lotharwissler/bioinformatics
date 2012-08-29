#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself
import fasta

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file" )
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
    if key == '-f':	args['fastafile'] = value
        
  if not args.has_key('fastafile'):
    stderr( "fasta file missing." )
    show_help()
  if not file_exists( args.get('fastafile') ):
    stderr( "fasta file does not exist." )
    show_help()
    
  return args

# =============================================================================
def get_sequences(file):
  seqcount, alnlength = 0, 0
  text = ''
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"):
      id = line[1:]
      if id.count(" ") > 0: id = id[:id.index(" ")]
      text += "\n" + id + "\n"
      seqcount += 1
    else:
      text += line
      if seqcount == 1: alnlength += len(line)
  fo.close()
  return text, seqcount, alnlength

# =============================================================================
# =============================================================================
def main( args ):
    for gid, seq in fasta.get_sequence_hash(args['fastafile']).iteritems():
      print string.join([gid, seq], "\t")

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
