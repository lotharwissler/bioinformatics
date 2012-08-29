#!/usr/bin/python
import os, sys
import string
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from tempfile import mkstemp
from collections import defaultdict



# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -a <path> -t <path> -m <N> -n <namespaces>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -a        annotation file in blast2go annot format (geneid <tab> goid)" )
  stdout( " -g        gene ontology obo file" )
  stdout( " -s        go slim obo" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "ha:g:s:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  map2slim = os.system("which map2slim > /dev/null")
  if not map2slim == 0: sys.exit("map2slim program not installed or in path. try 'sudo cpan GO::Parser'")

  args = {}
  for key, value in keys:
    if key == '-a': args['annot'] = value
    if key == '-g': args['goobo'] = value
    if key == '-s': args['slimobo'] = value
    
  if not args.has_key('annot'):
    stderr( "annot file argument missing." )
    show_help()
  elif not file_exists( args.get('annot') ):
    stderr( "annot file does not exist." )
    show_help()
  
  if not args.has_key('goobo'):
    stderr( "go obo file argument missing." )
    show_help()
  elif not file_exists( args.get('goobo') ):
    stderr( "go obo file does not exist." )
    show_help()

  if not args.has_key('slimobo'):
    stderr( "goslim obo file argument missing." )
    show_help()
  elif not file_exists( args.get('slimobo') ):
    stderr( "goslim obo file does not exist." )
    show_help()

  return args

# =============================================================================
def annot2fb(annot):
  gohash = {}
  fo = open(annot)
  fb = mkstemp('.fb', 'go2slim', '/tmp')[1]
  fw = open(fb, 'w')
  for line in fo:
    col = line.rstrip().split("\t")
    goid = col[1]
    if gohash.has_key(goid): continue
    gohash[goid] = 1
    fw.write(string.join([goid]*5 + [""]*12, "\t") + "\n")
  fo.close()
  fw.close()
  return fb

# =============================================================================
def annot2slim(annot, result):
  go2slim = defaultdict(list)
  fo = open(result)
  for line in fo:
    col = line.rstrip().split("\t")
    goid, slimid = col[0], col[4]
    go2slim[goid].append(slimid)
  fo.close()
  fo = open(annot)
  for line in fo:
    col = line.rstrip().split("\t")
    for slimid in go2slim[col[1]]:
      col[1] = slimid
      print string.join(col, "\t")
  fo.close()

# =============================================================================
def main(args):
  fb = annot2fb(args['annot'])
  result = mkstemp('.out', 'go2slim', '/tmp')[1]
  os.system("map2slim %s %s %s > %s" %(args['slimobo'], args['goobo'], fb, result))
  annot2slim(args['annot'], result)
  
  
# =============================================================================
args = handle_arguments()
main( args )
