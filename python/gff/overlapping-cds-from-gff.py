#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import tempfile
from low import *  # custom functions, written by myself
BEDBIN = "~/bin/bedtools/BEDTools-Version-2.14.3/bin/intersectBed"

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        gene feature file" )
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
    if key == '-f': args['gffile'] = value
    
  for key in ['gffile']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key):
      print >> sys.stderr, "missing argument", key
      show_help()
  return args

# =============================================================================
def statusbar(current, total, message="", width=40):
  progress = 1.0*current/total
  if message != "": message = "[" + message + "]"
  progressbar = "=" * int(progress*width)
  while len(progressbar) < width: progressbar += " " 
  sys.stderr.write("\r   0% " + progressbar + " 100% " + message)
  if progress == 1.0: sys.stderr.write("\n")
  

# =============================================================================
def find_overlapping_cds(gffile, overlap=20):
  fh, fn = tempfile.mkstemp()
  os.close(fh)
  cdsff, outfile = fn + ".gff", fn + ".out"
  if gffile.endswith(".gz"): readbin = "zcat"
  else: readbin = "cat"
  os.system("%s %s | awk -F \"\t\" '$3 == \"cds\" || $3 == \"CDS\" {print $0}' > %s" %(readbin, gffile, cdsff))
  os.system("cat %s | awk -F \"\t\" '$4 == 0 {$4 = 1}; {print}' > %s" %(cdsff, cdsff+'.1'))
  os.system("cat %s | awk -F \"\t\" '$4 < $5 {print $0}' > %s" %(cdsff+'.1', cdsff))
  os.system("cat %s | awk -F \"\t\" '$4 > $5 {print $1, $2, $3, $5, $4, $6, $7, $8, $9}' >> %s" %(cdsff+'.1', cdsff))
  os.system("%s -a %s -b %s -wo | awk -F \"\t\" '$7 != $16 && $19 >= %s {print $0}' > %s" %(BEDBIN, cdsff, cdsff, overlap, outfile))
  os.unlink(fn)
  os.unlink(cdsff)
  os.unlink(cdsff+'.1')
  return outfile

# =============================================================================
def find_parent(s):
  s = s[s.index("Parent=")+7:]
  s = s[:s.index(";")]
  return s

# =============================================================================
def parse_overlapping_gene_pairs(outfile):
  fo = open(outfile)
  for line in fo:
    cols = line.strip().split("\t")
    a1, a2 = cols[8], cols[17]
    p1, p2 = find_parent(a1), find_parent(a2)
    out = [p1, p2]
    out.sort()
    print string.join(out, "\t")
  os.unlink(outfile)

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  outfile = find_overlapping_cds(args['gffile'])
  parse_overlapping_gene_pairs(outfile)

# =============================================================================
args = handle_arguments()
main( args )

