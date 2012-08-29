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
	stdout( " -f        kegg html file" )
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
  statics = {}
  statics['entry'] = '^#ENTRY\s+(\S+)'
  statics['name'] = '^#NAME\s+(\S+)'
  statics['definition'] = '^#DEFINITION\s+(.*)$'
  oldlevel = ""
  hier = []
  for line in fo:
    for name, regex in statics.iteritems():
      if re.search( regex, line ):
        print "#%s\t%s" %(name, re.search( regex, line).group(1))

    if re.match( '[A-Z]\s+', line ):
      currentlevel = line[0]
      #print currentlevel
      rest = re.match( '[A-Z]\s+(.*)$', line ).group(1).strip()
      if not re.search( '\S+', rest ): continue
      rest = re.match( '(\S+)', rest ).group(1)
      if currentlevel > oldlevel:
        hier.append( strip_tags(rest) )
      elif currentlevel == oldlevel: 
        print string.join( hier, '/' )
        hier.pop()
        hier.append( strip_tags(rest) )
      else:
        hier.pop()
      
      oldlevel = currentlevel

  fo.close()


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
