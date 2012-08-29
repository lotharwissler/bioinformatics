#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself
import anydbm					# index databases (file hash)
from Bio import SeqIO # biopython stuff, to parse fasta files for instance

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -c <path> -o <path>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        blast.out file" )
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
		if key == '-f':	args['file'] = value
				
	if not args.has_key('file'):
		stderr( "blast.out file missing." )
		show_help()
	if not file_exists( args.get('file') ):
		stderr( "blast.out file does not exist." )
		show_help()
		
	return args



# =============================================================================
def parse_descr( text ):
  hash = {}
  if not re.search("GO:\d+.*evidence", text): 
    sys.stderr.write("return None.\n")
    return hash
  for match in re.finditer( '(GO:\d+)\s*\"([^"]+)\"\s*evidence', text ):
    id = match.group(1)
    description = match.group(2)
    hash[ id ] = description
  return hash
  

# =============================================================================
# =============================================================================
def main( args ):
  fo = open( args.get('file') )
  descr_index = None
  for line in fo:
    line = line.rstrip()
    cols = line.split("\t")
    if descr_index == None:
      for index, col in enumerate(cols):
        if re.search("GO:\d+", col):
          descr_index = index
          break
    descr = cols[ descr_index ]
    go_hash = parse_descr( descr )
    for goterm, godescr in go_hash.iteritems():
      L = []
      for index, col in enumerate(cols):
        if index == descr_index:
          L.append(goterm)
          L.append(godescr)
        else:
          L.append(col)
      print string.join(L,"\t")
  fo.close()
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
