#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself
from misa import MisaSSRspecies
import pickle

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print >> sys.stderr, "usage: " + sys.argv[0] + " -d <gff-folder>"
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -d        folder with gff files to parse" )
  stdout( " -f        misa file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hd:f:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-d': args['dir'] = value
    if key == '-f': args['misa'] = value
    
  if not args.has_key('dir'):
    print >> sys.stderr, "gff dir argument missing."
    show_help()
  elif not dir_exists( args.get('dir') ):
    print >> sys.stderr, "gff dir does not exist."
    show_help()

  if not args.has_key('misa'):
    print >> sys.stderr, "misa file argument missing."
    show_help()
  elif not file_exists( args.get('misa') ):
    print >> sys.stderr, "misa file does not exist."
    show_help()


  if not args['dir'].endswith("/"): args['dir'] += '/'
  return args

# =============================================================================
def get_ssrs(file):
    hash = defaultdict(list)
    fo = open(file)
    for line in fo: 
      if line.startswith("ID\t"): continue
      m = MisaSSRspecies(line)
      hash[m.species + '|' + m.geneid].append(m)
    fo.close()
    return hash

# =============================================================================
def get_features(dir):
  storage = ".misa.gff.storage.tmp"
  features = defaultdict(list)
  if not file_exists(storage):
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
      for line in fo: 
        if line.startswith("#") or len(line.rstrip()) == 0: continue
        columns = line.rstrip().split("\t")
        if len(columns) != 9: continue
        type = columns[2]
        if type != "gene" and type != "exon" and type != "intron": continue
        chr, start, stop, strand, descr = columns[0], columns[3], columns[4], columns[6], columns[8]
        key = string.join([species, chr], "|")
        features[key].append([type, int(start), int(stop)])
      fo.close()
      if gzip: os.system("gzip " + filename)

    fw = open(storage, "w")
    for key, features in features.iteritems():
      for feat in features:
        fw.write(string.join([key, feat[0], str(feat[1]), str(feat[2])], "\t") + "\n")
    fw.close()

  else:
    fo = open(storage)
    for line in fo:
      columns = line.rstrip().split("\t")
      key = columns[0]
      feat = list(columns[1:4])
      features[key].append(feat)
    fo.close()
  return features

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  SSRs = get_ssrs(args['misa']) 
  print >> sys.stderr, "misa output loaded..."
  Features = get_features(args['dir']) 
  print >> sys.stderr, "gff loaded..."
  total = len(SSRs)
  count = 0
  for key, ssrs in SSRs.iteritems():
    sys.stderr.write("\r   STATUS: %s/%s (%.2f%%)   CURRENT LOCATION: %s (%s SSRs, %s features)" %(count, total, (100.0*count/total), key, len(ssrs), len(Features[key])))
    for ssr in ssrs:
      type = 0
      for feat in Features[key]:
        ftype, fstart, fstop = feat[0], int(feat[1]), int(feat[2])
        if ssr.startpos >= fstart and ssr.startpos <= fstop:
          if not type or type == "gene": type = ftype
      if type == "gene": type = "UTR"
      if type == 0: type = "intergenic"
      print ssr.to_s() + "\t" + type
    count += 1
  sys.stderr.write("\n")



# =============================================================================
args = handle_arguments()
main( args )
