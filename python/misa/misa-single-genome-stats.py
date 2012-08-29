#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from misa import MisaSSRspecies


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        scaffold fasta file on which SSR identification was done" )
  stdout( " -m        misa out file (with species name as 1st column and localization feature in last column)" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:m:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-m': args['misa'] = value
    if key == '-f': args['fasta'] = value

  if not args.has_key('misa'):
    stderr( "misa file argument missing." )
    show_help()
  elif not file_exists( args.get('misa') ):
    stderr( "misa file does not exist." )
    show_help()
     
  if not args.has_key('fasta'):
    stderr( "fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('fasta') ):
    stderr( "fasta file does not exist." )
    show_help()
  
  return args

  
# =============================================================================
def get_scaffold_lengths(file):
  lengths = {}
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"):
      id = line[1:]
      if id.count(" ") > 0: id = id[:id.index(" ")]
      lengths[id] = 0
    else: lengths[id] += len(line.replace(" ", ''))
  return lengths

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  scaffold2length = get_scaffold_lengths(args['fasta'])
  scaffoldhash = {}
  fo = open(args['misa'])
  for line in fo:
    m = MisaSSRspecies(line)
    scaffold = m.geneid
    if not scaffoldhash.has_key(scaffold):
      scaffoldhash[scaffold] = {'ssrs':0, 'length':scaffold2length[scaffold], 'ssrlength':0, 'exonic':0, 'intronic':0, 'intergenic':0, "5'UTR":0, "3'UTR":0}
    scaffoldhash[scaffold]['ssrs'] += 1
    scaffoldhash[scaffold]['ssrlength'] += m.length
    if m.feature == 'E': scaffoldhash[scaffold]['exonic'] += 1
    elif m.feature == 'I': scaffoldhash[scaffold]['intronic'] += 1
    elif m.feature == '3': scaffoldhash[scaffold]["3'UTR"] += 1
    elif m.feature == '5': scaffoldhash[scaffold]["5'UTR"] += 1
    else: scaffoldhash[scaffold]['intergenic'] += 1
  fo.close()

  print string.join(["Scaffold", "Length (kb)", "SSRs (bp)", "# SSRs", "SSR.coverage", "SSRs/kb", "exonic", "intronic", "3'UTR", "5'UTR", "intergenic"], "\t")
  for scaffold, length in scaffold2length.iteritems():
    if not scaffoldhash.has_key(scaffold):
      sys.stdout.write(scaffold + "\t" + str(1.0*length/1000) + "\t")
      sys.stdout.write(string.join(["0"]*9, "\t") + "\n")
    else:
      hash = scaffoldhash[scaffold]
      sys.stdout.write(scaffold + "\t" + str(1.0*hash['length']/1000))
      sys.stdout.write("\t" + str(hash['ssrlength']))
      sys.stdout.write("\t" + str(hash['ssrs']))
      sys.stdout.write("\t" + "%0.2f%%" %(100.0*hash['ssrlength'] / hash['length']))
      sys.stdout.write("\t" + "%0.2f" %(1000.0*hash['ssrs'] / hash['length']))
      sys.stdout.write("\t" + str(hash['exonic']) + "\t" + str(hash['intronic']))
      sys.stdout.write("\t" + str(hash["3'UTR"]) + "\t" + str(hash["5'UTR"]))
      sys.stdout.write("\t" + str(hash['intergenic']) + "\n")



# =============================================================================
args = handle_arguments()
main( args )

