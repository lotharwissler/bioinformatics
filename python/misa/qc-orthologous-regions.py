#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import hashlib
from low import *  # custom functions, written by myself
from misa import MisaSSR
import newick
from collections import defaultdict
import pickle


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> " )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        orthologous region map" )
  stdout( " -g        all.gff" )
  stdout( " -m        gene2transcript2protein.map" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:g:m:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-g': args['gff'] = value
    if key == '-m': args['map'] = value
    
  if not args.has_key('file'):
    stderr( "orth.map file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "orth.map file does not exist." )
    show_help()
  if not args.has_key('map'):
    stderr( "gene2transcript2protein map file argument missing." )
    show_help()
  elif not file_exists( args.get('map') ):
    stderr( "gene2transcript2protein map file does not exist." )
    show_help()

  if not args.has_key('gff'):
    stderr( "gff file argument missing." )
    show_help()
  elif not file_exists( args.get('gff') ):
    stderr( "gff file does not exist." )
    show_help()


  return args


# =============================================================================
def coordinates_to_gene(file):
  hash = {}
  fo = open(file)
  for line in fo:
    cols = line.rstrip().split("\t")
    if not cols[3] == "gene": continue
    key = string.join([cols[0], cols[1], cols[4]], "|")
    value = [re.search("ID=([^;]+);", cols[9]).group(1), cols[7]]
    hash[key] = value
  fo.close()
  return hash


def gene_to_transcript(file):
  hash = {}
  fo = open(file)
  for line in fo:
    gid, tid = line.rstrip().split("\t")[0:2]
    hash[gid] = tid
  fo.close()
  return hash



def get_gene_features(file):
  exons = defaultdict(list)
  introns = defaultdict(list)
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    cols = line.split("\t")
    if cols[3] != "exon" and cols[3] != "intron": continue
    tid = re.search("Parent=([^;]+)", cols[9]).group(1)
    start, stop = cols[4], cols[5]
    strand = cols[7]
    if cols[3] == "exon": exons[tid].append([start, stop, strand])
    else: introns[tid].append([start, stop, strand])
#    hash[cols[0] + "|" + cols[1] + "|" + cols[3]] = 
  return exons, introns

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  gene2transcript = gene_to_transcript(args['map'])
  print >> sys.stderr, "gene2transcript loaded."
  coord2gene = coordinates_to_gene(args['gff'])
  print >> sys.stderr, "coord2gene loaded."
  exons, introns = get_gene_features(args['gff'])
  print >> sys.stderr, "gene features loaded."

  fo = open(args['file'])
  for line in fo:
    if line.startswith("#"): continue
    if len(line.split("\t")) < 9: continue
    line = line.rstrip()
    cols = line.split("\t")
    species1, species2 = cols[0:2]
    type = cols[2]
    chr1, chr2 = cols[3], cols[6]
    start1, start2 = cols[4], cols[7]
    stop1, stop2 = cols[5], cols[8]

    # remove regions with length=0 or where one region is significantly longer (150%)
    l1 = int(cols[5]) - int(cols[4])
    l2 = int(cols[8]) - int(cols[7])
    if l1 == 0 or l2 == 0: continue
    if float(max([l1,l2])) / float(min([l1,l2])) > 1.5 or (max([l1,l2]) - min([l1,l2])) > 5000: continue

    if type == "gene":
      key = string.join([species1, chr1, start1], "|")
      gid, strand1 = coord2gene[key]
      if not gene2transcript.has_key(gid): continue
      tid1 = gene2transcript[gid]
      exons1 = exons[tid1]
      introns1 = introns[tid1]

      key = string.join([species2, chr2, start2], "|")
      gid, strand2 = coord2gene[key]
      if not gene2transcript.has_key(gid): continue
      tid2 = gene2transcript[gid]
      exons2 = exons[tid2]
      introns2 = introns[tid2]

      if len(exons1) != len(exons2): continue
      if len(introns1) != len(introns2): continue

      cols.insert(6, strand1)
      cols.insert(10, strand2)

      # replace a gene by all its exons and introns
      for i in range(len(exons1)):
        ex1, ex2 = exons1[i], exons2[i]
        cols[2] = "exon"
        cols[4:7] = ex1
        cols[8:11] = ex2
        print string.join(cols, "\t")
      for i in range(len(introns1)):
        in1, in2 = introns1[i], introns2[i]
        cols[2] = "intron"
        cols[4:7] = in1
        cols[8:11] = in2
        print string.join(cols, "\t")
      continue
    
    key1 = string.join([species1, chr1, str(int(stop1) +1)], "|")
    key2 = string.join([species2, chr2, str(int(stop2) +1)], "|")
    gid, strand1 = coord2gene[key1]
    gid, strand2 = coord2gene[key2]
    cols.insert(6, strand1)
    cols.insert(10, strand2)
    print string.join(cols, "\t")

  fo.close()
#  print "exons equal:", ee, "exons unequal:", ue, "introns equal:", ei, "introns unqual:", ui, "both equal:", be

# =============================================================================
args = handle_arguments()
main( args )

