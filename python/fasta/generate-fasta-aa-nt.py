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
	stdout( "usage: " + sys.argv[0] + " -i <path> [-d <path>]" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -i        ID file" )
	stdout( " -d        directory to search for orthologs" )
	stdout( " " )
	
	sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hi:d:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-d':	args['dir'] = value
    if key == '-i':	args['idfile'] = value
        
  if args.has_key('dir') and not dir_exists( args.get('dir') ):
    stderr( "dir folder does not exist." )
    show_help()
  if not args.has_key('dir'): args['dir'] = './'
  if not args.get('dir').endswith('/'): args['dir'] = args.get('dir') + '/'
    
  if not args.has_key('idfile'):
    stderr( "id file missing." )
    show_help()
  if not file_exists( args.get('idfile') ):
    stderr( "id file does not exist" )
    show_help

  return args

  
# =============================================================================
# =============================================================================
def main( args ):
  idlist = read_from_file( args.get('idfile') ).splitlines()
  dir = args.get('dir')

  hash = {}
  for id in idlist:
    popenout = os.popen("grep -l \"%s\" %s*" %(id, dir))
    out = popenout.read()
    popenout.close()
    outlines = out.splitlines()

    hash[ id ] = outlines

  aafile = args.get('idfile') + '.aa'
  ntfile = args.get('idfile') + '.nt'
  for id,files in hash.iteritems():
    for file in files:
      if not file.endswith('.aa') and not file.endswith('.nt'): continue
      popenout = os.popen("grep -A 100 \"%s\" %s" %(id, file))
      out = popenout.read()
      popenout.close()
      outlines = out.splitlines()
      outlines.pop(0)
      
      if file.endswith('.aa'): outfile = aafile
      else: outfile = ntfile
      
      os.system( "echo \">%s\" >> %s" %( id, outfile ) )
      for line in outlines:
        if not line.startswith(">"): os.system( "echo \"%s\" >> %s" %( line, outfile ) )
        else: break

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
