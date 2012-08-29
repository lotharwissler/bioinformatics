#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print >> sys.stderr, "usage: " + sys.argv[0] + " -d <gff-folder>"
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -d        folder with gff files to parse" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hd:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-d': args['dir'] = value
    
  if not args.has_key('dir'):
    print >> sys.stderr, "gff dir argument missing."
    show_help()
  elif not dir_exists( args.get('dir') ):
    print >> sys.stderr, "gff dir does not exist."
    show_help()

  if not args['dir'].endswith("/"): args['dir'] += '/'
  return args


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  def process_gff_line(line, species):
    if line.startswith("#") or len(line.rstrip()) == 0: return
    columns = line.rstrip().split("\t")
    if len(columns) != 9: return
    type = columns[2]
    if type != "gene": return
    chr, start, stop, strand, descr = columns[0], columns[3], columns[4], columns[6], columns[8]
    id = re.search("ID=([^;]+);", descr).group(1)
    sys.stdout.write(species + "\t" + id + "\t")
    print string.join([chr, start, stop, strand], "\t")

# =============================================================================

  for filename in os.listdir(args['dir']):
    gzip = 0
    if not filename.endswith(".gff") and not filename.endswith(".gff.gz"): continue
    species = filename[:filename.index("-")]
    filename = args['dir'] +  filename
    if filename.endswith(".gff.gz"): gzip = 1
    if gzip: 
      os.system("gunzip " + filename)
      filename = filename[:-3]
    fo = open(filename)
    for line in fo: process_gff_line(line, species)
    fo.close()
    if gzip: os.system("gzip " + filename)



# =============================================================================
args = handle_arguments()
main( args )

