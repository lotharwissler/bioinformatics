#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from tempfile import mkstemp
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> -i" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file to import" )
  stdout( " -m        mapping: searchstring tab replacestring, one per line" )
  stdout( " -i        do in file, otherwise replace output to STDOUT" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:im:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  args['withinfile'] = 0
  for key, value in keys:
    if key == '-f': args['infile'] = value
    if key == '-m': args['mapfile'] = value
    if key == '-i': args['withinfile'] = 1
    
  if not args.has_key('infile'):
    stderr( "input file argument missing." )
    show_help()
  elif not file_exists( args.get('infile') ):
    stderr( "input file does not exist." )
    show_help()

  if not args.has_key('mapfile'):
    stderr( "map file argument missing." )
    show_help()
  elif not file_exists( args.get('mapfile') ):
    stderr( "map file does not exist." )
    show_help()
     
  return args


def get_map(mapfile):
  hash = {}
  fo = open(mapfile)
  for line in fo:
    line = line.strip()
    search, replace = line.split("\t")[0:2]
    hash[search] = replace
  fo.close()
  return hash

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  mapHash = get_map(args.get('mapfile'))
  if args.get('withinfile'):
    tmpfile = mkstemp(".tmp", "sr")[1]
    outstream = open( tmpfile, "w" )
  else:
    outstream = sys.stdout

  sout, serr = catch_bash_cmd_output( "wc -l %s" % args.get('infile') )
  total = int( sout.split()[0] )
  count = 0

  fo = open( args.get('infile') )
  for line in fo:
    line = line.rstrip()
    for s,r in mapHash.iteritems(): 
      if not s in line: continue
      line = line.replace(s,r)
    print >> outstream, line
    count += 1
    progress = int(50.0 * count / total) * "#"
    progress += (50 - len(progress)) * " "
    info("       0% " + progress + " 100%     ")
  fo.close()
  if args.get('withinfile'):
    outstream.close()
    os.system("mv %s %s" %( tmpfile, args.get('infile') ))
  info("       0% " + progress + " 100%    \n")


# =============================================================================
args = handle_arguments()
main( args )

