#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -s <path> -k \"regex\" -v \"regex\"" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -s        stockholm file" )
	stdout( " -k        regular expression for the key" )
	stdout( " -v        regular expression for the value" )
	stdout( " " )
	sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hs:k:v:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {}
	for key, value in keys:
		if key == '-s': args['stockholm'] = value
		if key == '-k':	args['keyregex'] = re.compile(value + '(.*)$' )
		if key == '-v':	args['valueregex'] = re.compile(value + '(.*)$' )
				
	if not args.has_key('keyregex'):
		stderr( "key regex missing." )
		show_help()
		
	if not args.has_key('valueregex'):
		stderr( "value regex missing." )
		show_help()

	if not args.has_key('stockholm'):
		stderr( "stockholm file missing." )
		show_help()
	if not file_exists( args.get('stockholm') ):
		stderr( "stockholm file does not exist." )
		show_help()
		
	return args

	
# =============================================================================
# =============================================================================
def main( args ):

	fo = open( args.get('stockholm') )
	kre = args.get('keyregex')
	vre = args.get('valueregex')
	key, value = '', ''
	for line in fo:
		if re.search( kre, line ):
			if key != '' and value != '':
				print ">%s" % key
				print value
				key, value = '', ''
			key = re.search( kre, line ).group(1).strip()
		if re.search( vre, line ):
			value = re.search( vre, line ).group(1).strip()
	fo.close()
	if key != '' and value != '':
		print ">%s" % key
		print value
	
	
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )