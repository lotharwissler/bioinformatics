#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import sqlite3
import glob
from low import *  # custom functions, written by myself


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        dir with fasta files on which SSR identification was done (*all-chromosome*.fasta)" )
  stdout( " -d        sqlite3 database" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:d:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-d': args['db'] = value
    if key == '-f': args['fasta'] = value

  if not args.has_key('db'):
    stderr( "db file argument missing." )
    show_help()
  elif not file_exists( args.get('db') ):
    stderr( "db file does not exist." )
    show_help()
     
  if not args.has_key('fasta'):
    stderr( "fasta dir argument missing." )
    show_help()
  elif not dir_exists( args.get('fasta') ):
    stderr( "fasta dir does not exist." )
    show_help()
  
  return args

  
# =============================================================================
def get_fasta_length(file):
  length = 0
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"): continue
    length += len(line.replace(" ", ''))
  return length


# =============================================================================
def get_length_hash(fastadir, species):
  print >> sys.stderr, "getting genome sizes..."
  fastahash = {}
  for spec in species:
    dest = fastadir + '/' + spec + '*all-chromosome*.fasta'
    files = glob.glob(dest)
    if len(files) == 0: 
      sys.exit("no fasta file found for species %s (%s). aborting." %(spec, dest))
    elif len(files) > 1:
      sys.exit("more than 1 file found for species %s (%s). aborting." %(spec, dest))
    else:
      fastahash[spec] = get_fasta_length(files[0])
  return fastahash
  

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  # TODO: swtich reading to the sqlite3 db
  # execute appropriate queries including the gene_features table
  conn = sqlite3.connect(args['db'])
  species = [str(e[0]) for e in conn.execute("SELECT DISTINCT species FROM ssrs").fetchall()]
  lengthhash = get_length_hash(args['fasta'], species)

  sys.stdout.write(string.join(["Species", "Length (bp)", "SSRs (bp)", "# SSRs", "SSR.coverage", "SSRs/kb"], "\t") + "\t")
  sys.stdout.write(string.join(["genomic." + str(e) for e in range(1,7)], "\t") + "\t")
  sys.stdout.write(string.join(["exonic." + str(e) for e in range(1,7)], "\t") + "\t")
  sys.stdout.write(string.join(["intronic." + str(e) for e in range(1,7)], "\t") + "\t")
  sys.stdout.write(string.join(["intergenic." + str(e) for e in range(1,7)], "\t") + "\t")
  sys.stdout.write("\n")

  for spec in species:
    c = conn.execute("SELECT length FROM ssrs WHERE species='%s'" % spec)
    ssrlength = 0
    rows = c.fetchmany()
    while rows:
      for l in [e[0] for e in rows]: ssrlength += l
      rows = c.fetchmany()
    ssrcount = conn.execute("SELECT COUNT(id) FROM ssrs WHERE species='%s'" % spec).fetchall()[0][0]
    sys.stdout.write(spec + "\t")
    sys.stdout.write(str(lengthhash[spec]) + "\t")
    sys.stdout.write(str(ssrlength) + "\t")
    sys.stdout.write(str(ssrcount) + "\t")
    sys.stdout.write(str(1.0*ssrlength/lengthhash[spec]) + "\t")
    sys.stdout.write(str(1.0*ssrcount/lengthhash[spec]*1000) + "\t")
    ssrtypefreq = {}
    for row in conn.execute("SELECT LENGTH(motif),COUNT(id) FROM ssrs WHERE species='%s' GROUP BY LENGTH(motif)" % spec).fetchall(): ssrtypefreq[str(row[0])] = row[1]
    for i in range(1,7): sys.stdout.write(str(ssrtypefreq[str(i)]) + "\t")
    ssrtypefreq = {}
    for row in conn.execute("SELECT LENGTH(motif),COUNT(id) FROM ssrs WHERE species='%s' AND gene_feature='exon' GROUP BY LENGTH(motif)" % spec).fetchall(): ssrtypefreq[str(row[0])] = row[1]
    for i in range(1,7): sys.stdout.write(str(ssrtypefreq[str(i)]) + "\t")
    ssrtypefreq = {}
    for row in conn.execute("SELECT LENGTH(motif),COUNT(id) FROM ssrs WHERE species='%s' AND gene_feature='intron' GROUP BY LENGTH(motif)" % spec).fetchall(): ssrtypefreq[str(row[0])] = row[1]
    for i in range(1,7): sys.stdout.write(str(ssrtypefreq[str(i)]) + "\t")
    ssrtypefreq = {}
    for row in conn.execute("SELECT LENGTH(motif),COUNT(id) FROM ssrs WHERE species='%s' AND gene_feature IS NULL GROUP BY LENGTH(motif)" % spec).fetchall(): ssrtypefreq[str(row[0])] = row[1]
    for i in range(1,7): sys.stdout.write(str(ssrtypefreq[str(i)]) + "\t")
    sys.stdout.write("\n")

  conn.close()

# =============================================================================
args = handle_arguments()
main( args )

