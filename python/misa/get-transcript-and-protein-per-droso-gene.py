#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself
from misa import MisaSSR
import pickle

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print >> sys.stderr, "usage: " + sys.argv[0] + " -f <fasta>"
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        all.translations.fasta" )
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
    if key == '-f': args['fasta'] = value
    
  if not args.has_key('fasta'):
    print >> sys.stderr, "fasta file argument missing."
    show_help()
  elif not file_exists( args.get('fasta') ):
    print >> sys.stderr, "fasta file does not exist."
    show_help()


  return args

# =============================================================================
def get_ssrs(file):
    hash = defaultdict(list)
    fo = open(file)
    for line in fo: 
      if line.startswith("ID\t"): continue
      m = MisaSSR(line)
      hash[m.geneid].append(m)
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

  geneHash = {}
  fo = open(args['fasta'])
  for line in fo:
    if not line.startswith(">"): continue
    proteinid = re.search(" ID=(\S+);", line).group(1)
    geneid, transcriptid = re.search(" parent=(\S+);", line).group(1).split(",")[0:2]
    length = int(re.search(" length=(\d+)", line).group(1))
    if not geneHash.has_key(geneid) or geneHash[geneid]['length'] < length: 
      geneHash[geneid] = {'protein':proteinid, 'transcript':transcriptid, 'length':length}
  fo.close()

  for geneid, hash in geneHash.iteritems():
    print string.join([geneid, hash['transcript'], hash['protein']], "\t")



# =============================================================================
args = handle_arguments()
main( args )
