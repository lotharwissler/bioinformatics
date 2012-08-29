#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import math        # match functions
from low import *  # custom functions, written by myself

REGEX = re.compile("(\d+)$")

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        text flat file to analyze" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
    
  return args


# =============================================================================
def is1higherthan( text1, text2, regex=REGEX ):

  def splittext( text, regex ):
    return regex.split( text )[0], int(regex.split( text )[1])

  id1, number1 = splittext( text1, regex )
  id2, number2 = splittext( text2, regex )
  if id1 != id2: return 0 
  if (number1 +1) == number2: return 1
  return 0

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  fo = open( args.get('file') )
  lines = fo.readlines()
  fo.close()

  started_at = ""

  for i in range(1,len(lines)):
    line0, line1 = lines[i-1], lines[i]
    if started_at == "": started_at = line0
    if i < (len(lines)-1) and is1higherthan( line0, line1 ): continue
    print string.join([started_at.rstrip(), line0.rstrip()], "\t")
    started_at = ""
    
# =============================================================================
args = handle_arguments()
main( args )

