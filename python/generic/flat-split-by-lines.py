#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from low import *

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print "splits a flat file into chunks. Options: (1) N number of lines per chunk. (2) N number of chunks of equal size"
  print "usage: " + sys.argv[0] + " -f <file> [-i <chunks> -l <lines>]"
  print " "
  print " option    description"
  print " -h        help (this text here)"
  print " -f        flat file to split"
  print " -l        number of lines per chunk"
  print " -i        number of equally sized chunks"
  print " "
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    print >> sys.stderr, "no arguments provided."
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:i:l:" )
  except getopt.GetoptError:
    print >> sys.stderr, "invalid arguments provided."
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-l': args['l'] = int(value)
    if key == '-i': args['i'] = int(value)
    
  if not args.has_key('file'):
    print >> sys.stderr, "import file argument missing."
    show_help()
  elif not file_exists( args.get('file') ):
    print >> sys.stderr, "import file does not exist."
    show_help()

  if not args.has_key('l') and not args.has_key('i'):
    print >> sys.stderr, "l or i missing."
    show_help()
    
  return args


def get_number_of_lines(file):
  lines = 0
  fo = open(file)
  for line in fo: lines += 1
  return lines

# =============================================================================
def get_lines_in(ifile):
  lc = 0
  fo = open(ifile)
  for line in fo: lc += 1
  fo.close()
  return lc
  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  totallines = get_lines_in(args.get('file'))
  linecount, filecount = 0, 1
  if args.has_key('i'): rotate = int(math.ceil( 1.0 * totallines / args.get('i') ))
  else: rotate = args.get('l')

  digits = len(str(math.ceil(1.0*totallines/rotate)))
  fw = open( args.get('file') + '.' + add_leading_zeroes(filecount, digits), 'w' )
  fo = open( args.get('file') )
  for line in fo:
    linecount += 1
    if ((linecount % rotate) == 1 and linecount > 1) or (rotate == 1 and linecount > 1):
      filecount += 1
      fw.close()
      fw = open( args.get('file') + '.' + add_leading_zeroes(filecount, digits), 'w' )
    fw.write(line)    
  fo.close()
  fw.close()
  

# =============================================================================
args = handle_arguments()
main( args )

