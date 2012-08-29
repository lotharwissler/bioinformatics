#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
import gff3

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -i        surrounding genes (ID1,ID2)" )
  stdout( " -g        gff file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """

  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hi:g:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-i': args['genes'] = value.split(',')
    if key == '-g': args['gffile'] = value
    
  for key in ['genes', 'gffile']:
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
def get_coordinates_for_diagram(gfhash, genes):
  positions = []
  for scaffold, gfs in gfhash.iteritems():
    for gf in gfs:
      if gf.ftype != 'mRNA' or not gf.get_attributes().has_key('ID'): continue
      if gf.get_attributes()['ID'] in genes: 
        positions += [gf.start, gf.stop]
        print gf.start, gf.stop
      if len(positions) == 4: return scaffold, min(positions), max(positions)
        
    
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  gfhash = gff3.get_gff_hash(args['gffile'])
  gid, Diagstart, Diagstop = get_coordinates_for_diagram(gfhash, args['genes'])
  print gid, Diagstart, Diagstop
  
  gd_diagram = GenomeDiagram.Diagram(gid)
  gd_track_for_features = gd_diagram.new_track(1, name="Annotated Genes")
  gd_feature_set = gd_track_for_features.new_set()
  
  for gf in sorted(gfhash[gid], key=lambda x: x.start):
    if gf.ftype != 'mRNA': continue
    if gf.start > gf.stop: gf.start, gf.stop = gf.stop, gf.start
    if gf.stop < Diagstart or gf.start > Diagstop: continue
    f = SeqFeature(FeatureLocation(max([gf.start, Diagstart]), min([gf.stop, Diagstop])), strand=int(gf.strand+'1'), type=gf.get_attributes()['ID'])
    gd_feature_set.add_feature(f, label=True, label_size=10, label_angle=0, sigil="ARROW")
    print gf.get_attributes()['ID'], gf.start, gf.stop
  
  gd_diagram.draw(start=Diagstart, end=Diagstop, format='linear', fragments=1, pagesize=(100*cm, 4*cm))
  outfile = os.path.split(args['gffile'])[1] + "_"+ string.join(args['genes'], "_") + '.pdf'
  gd_diagram.write(outfile, "PDF")
  print outfile

# =============================================================================
args = handle_arguments()
main( args )

