#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from collections import defaultdict
from xml.dom import minidom


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        go term obo-xml file" )
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
    if key == '-f': args['obo'] = value
    
  if not args.has_key('obo'):
    stderr( "obo file argument missing." )
    show_help()
  elif not file_exists( args.get('obo') ):
    stderr( "obo file does not exist." )
    show_help()
  
  return args

# =============================================================================
class GOTerm():

  def __init__(self, xml):
    self.id = xml.getElementsByTagName("id")[0].firstChild.data
    self.name = xml.getElementsByTagName("name")[0].firstChild.data
    self.namespace = xml.getElementsByTagName("namespace")[0].firstChild.data
    self.alt_ids = [node.firstChild.data for node in xml.getElementsByTagName("alt_id")]

# =============================================================================
def read_obo( file ):
  hash = {}
  xmldoc = minidom.parse(file)
  for term in xmldoc.getElementsByTagName('term'):
    goterm = GOTerm(term)
    hash[goterm.id] = goterm
    for alt_id in goterm.alt_ids: 
      if not hash.has_key(alt_id): hash[alt_id] = goterm
  print >> sys.stderr, "goterms read from obo: %s" % len(hash)
  return hash

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  gohash = read_obo(args['obo'])
  for goid, goterm in gohash.iteritems():
    print goid + "\t" + goterm.name

# =============================================================================
args = handle_arguments()
main( args )

