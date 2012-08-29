#!/usr/bin/python
import os, sys 				# low level handling, such as command line stuff
import getopt					# comand line argument handling
import anydbm					# index databases (file hash)
from low import *			# collection of generic self-defined functions


# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> -o <path>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        xdom file" )
	stdout( " -i        indexed ndb file" )
	stdout( " -n        column to look up [0..n]" )
	stdout( " " )
	sys.exit(1)
	
# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
  	stderr( "no arguments provided." )
  	show_help()	
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:i:n:" )
  except getopt.GetoptError:
  	stderr( "invalid arguments provided." )
  	show_help()
  
  args = {}
  for key, value in keys:
    if key == '-f': args['xdom'] = value
    if key == '-i':	args['dbm'] = value
    if key == '-n':	args['column'] = int(value)
  
  if not args.has_key('xdom'):
  	stderr( "xdom file missing." )
  	show_help()
  if not file_exists( args.get('xdom') ):
  	stderr( "xdom file does not exist." )
  	show_help()
  	
  if not args.has_key('dbm'):
  	stderr( "dbm file missing." )
  	show_help()
  if not file_exists( args.get('dbm') ):
  	stderr( "dbm file does not exist." )
  	show_help()
  
  if not args.has_key('column'):
    stderr( "column index missing." )
    show_help()
  
  return args

	
# =============================================================================
# =============================================================================
def main( args ):
  DBM = anydbm.open( args.get('dbm'), 'r' )
  fo = open( args.get('xdom') )
  n = args.get('column')
  key, value = '', ''
  for line in fo:
    line = line.rstrip()
    if line.endswith('\n'): line = line.replace('\n','')
    if line.startswith('>'):
      print line
  		#if key != '' and value != '':
  		#	sys.stdout.write( ">%s\n%s" %(key,value) )
  		#	key, value = '', ''
  		#key = line[1:].rstrip()
    else:
      value = line.rstrip()
      pid = value.split()[ n ]
      if not DBM.has_key( pid ):
      	print "DBM does not contain the following key:", pid
      else: value += "\t" + DBM.get(pid)
      print value 
  fo.close()
  #if key != '' and value != '':
  #	sys.stdout.write( ">%s\n%s" %(key,value) )
  DBM.close()

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
args = handle_arguments(  )
main( args )

