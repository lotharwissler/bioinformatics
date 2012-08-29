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
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        kegg KO annotation file" )
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
    if key == '-f':	args['file'] =  value
        
  if not args.has_key('file'):
    stderr( "kegg file missing." )
    show_help()
  if not file_exists( args.get('file') ):
    stderr( "kegg file does not exist." )
    show_help()

  return args


# =============================================================================
def strip_tags(value):
  "Return the given HTML with all tags (+ KEGG tags) stripped."
  value = re.sub(r'<[^>]*?>', '', value)
  value = re.sub(r'\[.*\]', '', value)
  return value


def read_KOs( file ):

  def next_entry(fo):
    pathlist = []
    definition = ""
    line = fo.readline().rstrip()
    if line == '': 
      return fo, None, None
    entry = re.match('^ENTRY\s+(\S+)', line).group(1)
    line = fo.readline().rstrip()
    line = fo.readline().rstrip()
    if re.match( '^DEFINITION\s+(.*)$',line):
      definition = re.search( '^DEFINITION\s+(.*)$', line ).group(1)
      line = fo.readline().rstrip()
    while line.startswith('CLASS') or line.startswith(' '):
      if re.search('\[\S+:\S+\]', line):
        pathlist.append( re.search('\[(\S+:\S+)\]',line).group(1) )
      line = fo.readline().rstrip()
      
    while line != '///':
      line = fo.readline().rstrip()

    if definition != "": entry += "\t" + definition
    return fo, entry, pathlist
  
  fo = open( file )
  kohash = {}
  while 1:
    fo, id, pathlist = next_entry( fo )
    if id == None: break
    print ">%s\n%s" %(id, string.join(pathlist,"\t"))
    
  fo.close()

# =============================================================================
# =============================================================================
def main( args ):
  
  kohash = read_KOs( args.get('file') )


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
