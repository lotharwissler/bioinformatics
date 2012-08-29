#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import math        # match functions
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> -i -n" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        tab delimited input file" )
  stdout( " -1        keep first column" )
  stdout( " -r        regex for the column header to mark as to keep" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:r:1" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  args['keepfirstcol'] = 0
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-r': args['regex'] = re.compile(value)
    if key == '-1': args['keepfirstcol'] = 1
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()

  if not args.has_key('regex'):
    stderr( "regex argument missing." )
    show_help()
  
  return args


# =============================================================================
def get_header( file ):
  fo = open(file)
  header = fo.readline().rstrip()
  fo.close()
  return header


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  headline = get_header( args.get('file') )
  columns = headline.split("\t")
  regex = args.get('regex')
  keepindices = []
  for i in range(len(columns)):
    if regex.search(columns[i]): 
      keepindices.append(i)
      #sys.stderr.write("marked:\t%02d\t%s\n" % (i, columns[i]))
    elif i == 0 and args.get('keepfirstcol'): 
      keepindices.append(i)
      #sys.stderr.write("marked:\t%02d\t%s\n" % (i, columns[i]))

  fo = open(args.get('file'))
  for line in fo:
    line = line.rstrip()
    columns = line.split("\t")
    out = []
    for i in keepindices: out.append( columns[i] )
    print string.join(out, "\t")
  fo.close()

# =============================================================================
args = handle_arguments()
main( args )

