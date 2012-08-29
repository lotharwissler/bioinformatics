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
	stdout( " -f        input fasta file" )
	stdout( " -o        output dbm file" )
	stdout( " " )
	sys.exit(1)
	
# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hf:o:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {}
	for key, value in keys:
		if key == '-f': args['fasta'] = value
		if key == '-o':	args['out'] = value
	
	if not args.has_key('fasta'):
		stderr( "fasta file missing." )
		show_help()
	if not file_exists( args.get('fasta') ):
		stderr( "fasta file does not exist." )
		show_help()
	
	if not args.has_key('out'):
		stderr( "out file missing." )
		show_help()
		
	return args

	
# =============================================================================
# =============================================================================
def main( args ):
  DBM = anydbm.open( args.get('out'), 'c' )
  sout, serr = catch_bash_cmd_output( "grep '>' -c %s" %args.get('fasta') )
  total = int( sout )
  added = 0
  fo = open( args.get('fasta') )
  key, value = '', ''
  for line in fo:
    line = line.rstrip()
    if line.startswith('>'):
      if key != '' and value != '':
        #print key + "\t" + value
        added += 1
        DBM[ key ] = value
        sys.stderr.write('\r\tindexing:\t%s\t%01.2f%%' %(added,100.0*added/total) )
        sys.stderr.flush()
        key, value = '', ''
      key = re.match(">(\S+)", line).group(1)
    else:
      value += line.rstrip()
  fo.close()
  if key != '' and value != '':
    added += 1
    DBM[ key ] = value
    #print key + "\t" + value
  DBM.close()
  sys.stderr.write('\r\tindexing:\t%s\t%01.2f%%\ndone.\n' %(added,100.0*added/total) )

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
args = handle_arguments(  )
main( args )

