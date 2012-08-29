#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
import math
from low import *			# custom functions, written by myself

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> -n <x> -i <x>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        fasta file" )
	stdout( " -n        size of each new fasta file (# seq)" )
	stdout( " -i        number of fasta files to split into" )
	stdout( " " )
	sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hf:n:i:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {}
	for key, value in keys:
		if key == '-f': args['fasta'] = value
		if key == '-n':	args['n'] = int(value)
		if key == '-i':	args['i'] = int(value)
				
	if not args.has_key('n') and not args.has_key('i'):
		stderr( "n or i missing." )
		show_help()

	if not args.has_key('fasta'):
		stderr( "fasta file missing." )
		show_help()
	if not file_exists( args.get('fasta') ):
		stderr( "fasta file does not exist." )
		show_help()
		
	return args

	
# =============================================================================
# =============================================================================
def main( args ):
  sout, serr = catch_bash_cmd_output( "grep '>' -c %s" % args.get('fasta') )
  total = int( sout )
  cut = total
  seqcount = 0
  filecount = 1

  if args.has_key('i'): cut = int(math.ceil( 1.0 * total / args.get('i') ))
  else: cut = args.get('n')


  fw = open( args.get('fasta') + '.' + add_leading_zeroes(filecount, 6), 'w' )
  handle = open(args.get('fasta'))
  for line in handle:

    if line[0] == ">": 
      seqcount += 1
      if ((seqcount % cut) == 1 and seqcount > 1) or (cut == 1 and seqcount > 1):
        filecount += 1
        fw.flush()
        fw.close()
        fw = open( args.get('fasta') + '.' + add_leading_zeroes(filecount, 6), 'w' )

    fw.write(line)

  fw.flush()
  fw.close()
  infomsg( "total.seq.count: %s | split.count: %s | file.count: %s" %(total, cut, filecount) )
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
