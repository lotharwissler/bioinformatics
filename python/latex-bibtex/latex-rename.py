#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "renames files so that they can be included in latex documents." )
  stdout( "this means that all dots are removed except for the last one of the actual file extension." )
  stdout( "dots are replaced by \"_\" by default." )
  stdout( "usage: " + sys.argv[0] + " -f <path> [-r <x>]" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file" )
  stdout( " -r        replace dot with this sign (default: \"_\")" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hf:r:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {}
	for key, value in keys:
		if key == '-f': args['file'] = value
		if key == '-r':	args['r'] = str(value)
				
	if not args.has_key('file'):
		stderr( "file missing." )
		show_help()
	if not file_exists( args.get('file') ):
		stderr( "file does not exist." )
		show_help()
		
	return args

	
# =============================================================================
# =============================================================================
def main( args ):
  oldfilename = args.get('file')
  path, filename = os.path.split(oldfilename)
  base, ext = os.path.splitext(filename)
  if args.has_key('r'): r = args.get('r')
  else: r = '_'
  base = base.replace('.',r)
  if path != "":
    newfilename = path + '/' + base + ext
  else:
    newfilename = base + ext
  os.system( "mv %s %s" %(oldfilename, newfilename) )
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
