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
  stdout( " -f        kegg html file" )
  stdout( " -k        kegg KO annotation file" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:k:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f':	args['file'] =  value
    if key == '-k':	args['ko'] = value
        
  if not args.has_key('file'):
    stderr( "kegg file missing." )
    show_help()
  if not file_exists( args.get('file') ):
    stderr( "kegg file does not exist." )
    show_help()

  if not args.has_key('ko'):
    stderr( "ko file missing." )
    show_help()
  if not file_exists( args.get('ko') ):
    stderr( "ko file does not exist." )
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
    hash = {}
    hash['class'] = ""
    line = fo.readline().rstrip()
    if line == '': 
      return fo, None, None
    entry = re.match('^ENTRY\s+(\S+)', line).group(1)
    line = fo.readline().rstrip()
    hash['name'] = re.match('^NAME\s+(.*)$',line).group(1)
    line = fo.readline().rstrip()
    if re.match( '^DEFINITION\s+(.*)$',line):
      hash['definition'] = re.match('^DEFINITION\s+(.*)$',line).group(1)
      line = fo.readline().rstrip()
    while line.startswith('CLASS') or line.startswith(' '):
      hash['class'] += line[12:]
      line = fo.readline().rstrip()
      
    while line != '///':
      line = fo.readline().rstrip()

    return fo, entry, hash
  
  fo = open( file )
  kohash = {}
  while 1:
    fo, id, hash = next_entry( fo )
    if id == None: break
    kohash[id] = hash
    
  fo.close()
  return kohash

# =============================================================================
# =============================================================================
def main( args ):
  
  kohash = read_KOs( args.get('ko') )
  regex = re.compile('[A-Z]\s+(\S+);\s+<a href=".*">(\S+)</a>')
  fo = open( args.get('file'), 'r' )
  for line in fo:
    if re.match( regex, line ):
      sid, koid = re.match( regex, line ).groups()
      print string.join( [sid, koid, kohash[koid].get('name'), kohash[koid].get('class')], "\t" )

  fo.close()


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
