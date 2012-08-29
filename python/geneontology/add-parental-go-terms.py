#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from goterm import GOTerm
from collections import defaultdict


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> -i -n" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -a        go annot file (geneid <tab> goid)" )
  stdout( " -o        go obo file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "ha:o:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-a': args['annot'] = value
    if key == '-o': args['obo'] = value
    
  if not args.has_key('annot'):
    stderr( "annot file argument missing." )
    show_help()
  elif not file_exists( args.get('annot') ):
    stderr( "annot file does not exist." )
    show_help()
  
  if not args.has_key('obo'):
    stderr( "obo file argument missing." )
    show_help()
  elif not file_exists( args.get('obo') ):
    stderr( "obo file does not exist." )
    show_help()
  
  return args

  
# =============================================================================
def read_obo( file ):
  hash = {}
  fo = open(file)
  while 1:
    line = fo.readline()
    if line.startswith("[Typedef]"): break
    if not line.startswith("[Term]"): continue
    obolines = []
    while 1:
      line = fo.readline()
      if line.count(":") == 0: break
      obolines.append(line)
    goterm = GOTerm(obolines)
    hash[goterm.get_id()] = goterm
    for alt_id in goterm.get_alt_ids():
      hash[alt_id] = goterm
  fo.close()
  print >> sys.stderr, "goterms read from obo: %s" % len(hash)
  return hash


def get_parents_of_goterms(file, gohash):

  def get_parents(goid):
    parents = gohash[goid].get_is_a_goids()
    return parents + [get_parents(p) for p in parents]

  def flatten(x):
    result = []
    for el in x:
      if hasattr(el, "__iter__") and not isinstance(el, basestring): result.extend(flatten(el))
      else: result.append(el)
    return result

  hash = {}
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    geneid, gotermid = line.split("\t")[0:2]
    if not gotermid.startswith("GO:"): continue
    if hash.has_key(gotermid): continue
    hash[gotermid] = []
  fo.close()
  allannotgotermids = hash.keys()
  for goid in allannotgotermids: hash[goid] = list(set(flatten(get_parents(goid))))
  return hash


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  print >> sys.stderr, "reading obo file ..."
  gohash = read_obo(args['obo'])
  print >> sys.stderr, "gather parents for all annotated goterms ..."
  goterm2parents = get_parents_of_goterms(args['annot'],gohash)

  print >> sys.stderr, "producing output ..."
  fo = open(args['annot'])
  for line in fo:
    line = line.rstrip()
    geneid, goid = line.split("\t")[0:2]
    if not goid.startswith("GO:"): continue
    print geneid + "\t" + goid
    for id in goterm2parents[goid]:
      print geneid + "\t" + id
  fo.close()
  print >> sys.stderr, "done."

# =============================================================================
args = handle_arguments()
main( args )

