#!/usr/bin/python
import os, sys 				# low level handling, such as command line stuff
import getopt					# comand line argument handling
import anydbm					# index databases (file hash)
from low import *			# collection of generic self-defined functions


# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> -e <string> [-i <n>]" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        table file" )
	stdout( " -i        table column index containing the lookup name [default: 0]" )
	stdout( " -c        annotation file column to use [default: all]" )
	stdout( " -l        annotation file line(s) to use [default: first]" )
	stdout( " -e        file extension to look for (= lookupname.extension)" )
	stdout( " " )
	sys.exit(1)
	
# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
  	stderr( "no arguments provided." )
  	show_help()	
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:i:e:" )
  except getopt.GetoptError:
  	stderr( "invalid arguments provided." )
  	show_help()
  
  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-i':	args['col'] = int( value )
    if key == '-e':	args['ext'] = value
  
  if not args.has_key('file'):
  	stderr( "table file missing." )
  	show_help()
  if not file_exists( args.get('file') ):
  	stderr( "table file does not exist." )
  	show_help()
  	
  if not args.has_key('col'):
     args['col'] = 0
  
  if not args.has_key('ext'):
    stderr( "file extension missing." )
    show_help()
  
  return args

	
# =============================================================================
# =============================================================================
def main( args ):
  fo = open( args.get('file') )
  for line in fo:
    line = line.rstrip()
    columns = line.split("\t")
    lookup = columns[ args.get('col') ]
    lookupfile = lookup + args['ext']
    if file_exists( lookupfile):
      ft = open( lookupfile )
      lines = ft.readlines()
      ft.close()
      # TODO:
      # get lines, get column
      # then add to table
      # print the new line

  fo.close()

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
args = handle_arguments(  )
main( args )

