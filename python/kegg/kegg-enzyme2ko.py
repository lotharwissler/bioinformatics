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
	stdout( "usage: " + sys.argv[0] + " -f <path>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        kegg ko file file" )
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
		stderr( "kegg file missing." )
		show_help()
	if not file_exists( args.get('file') ):
		stderr( "kegg file does not exist." )
		show_help()
		
	return args


# =============================================================================
def strip_tags(value):
  "Return the given HTML with all tags (+ KEGG tags) stripped."
  value = re.sub(r'<[^>]*?>', '', value)
  value = re.sub(r'\[.*\]', '', value)
  return value

# =============================================================================
# =============================================================================
def main( args ):
  fo = open( args.get('file'), 'r' )
  ko_regex = re.compile( "^ENTRY\s+(K\S+)" )
  enzyme_regex = re.compile( "\s+EC:\s+([0-9.]+)" )

  ko, enzyme = "", ""
  for line in fo:
    line = line.rstrip()
    if line.startswith("///"): 
      ko, enzyme = "", ""
      continue
    if ko == "":
      if re.search( ko_regex, line): ko = re.search( ko_regex, line ).group(1)
    else:
      if re.search( enzyme_regex, line):
        enzyme = re.search( enzyme_regex, line ).group(1)
        print "%s\t%s" % ( ko, enzyme )

  fo.close()


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
