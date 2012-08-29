#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from collections import defaultdict
import glob
import gff3


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> ..." )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        flat file" )
  stdout( " -s        separator of flat file (default: tab)" )
  stdout( " -a        action [INSERT|UPDATE]" )
  stdout( " -t        sql table name" )
  stdout( " " )
  stdout( " Field names are extracted from the header line (first line, must start with #)." )
  stdout( " NULL named fields are irgnored, the rest gets imported." )
  stdout( " UPDATES are only possible with a given ID. Thus, the header must contain " )
  stdout( " a column named ID which will be used to generate an UPDATE ... WHERE id='ID'." )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """

  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:a:s:t:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'separator':"\t", 'action':"INSERT"}
  for key, value in keys:
    if key == '-f': args['flatfile'] = value
    if key == '-s': args['separator'] = value
    if key == '-a': args['action'] = value.upper()
    if key == '-t': args['table'] = value
    
  for key in ['flatfile', 'separator', 'action', 'table']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key): show_help()
  return args

# =============================================================================
def get_blastout_hash(file):
  hash = defaultdict(int)
  fo = open(file)
  for line in fo:
    qid = line.split("\t")[0]
    hash[qid] = 1
  fo.close()
  return hash

# =============================================================================
def gather_blast_output(bdir):
  hash = {}
  for filename in glob.glob(bdir + '/*.blastout'):
    s = os.path.split(filename)[1][:4]
    hash[s] = get_blastout_hash(filename)
  return hash

# =============================================================================
def get_scaffolds(file):

  def add_feature(hash, gf):
    if not hash.has_key(gf.seqid): hash[gf.seqid] = {}
    hash[gf.seqid][gf.start] = gf.get_attributes()['ID']
    return hash

  hash = {}
  fo = open(file)
  for line in fo:
    if line.startswith("#"): continue
    gf = gff3.GeneFeature(line)
    if gf.ftype != "mRNA": continue
    hash = add_feature(hash, gf)
  fo.close()

  outhash = {}
  for scaffold, h in hash.iteritems():
    outhash[scaffold] = [h[key] for key in sorted(h.iterkeys())]

  return outhash

# =============================================================================
def gather_genes_on_scaffolds(gffdir):
  hash = {}
  for filename in glob.glob(gffdir + '/*.gff'):
    s = os.path.split(filename)[1][:4]
    hash[s] = get_scaffolds(filename)
  return hash

# =============================================================================
def get_neighbors(pid, geneids):
  index = geneids.index(pid)
  left = geneids[max([index-3,0]):index]
  right = geneids[index+1:min([index+4,len(geneids)])]
  return (left, right)

# =============================================================================
def escape4sql(string):
  if string.count("'") == 0: return string
  return string.replace("'", "\\'")

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  fo = open(args['flatfile'])
  action, table, fieldnames = args['action'], args['table'], ""
  print "SET autocommit=0;"
  for line in fo:
    if line.startswith("#"):
      if fieldnames == "": fieldnames = [e.strip().upper() for e in line[1:].split(args['separator'])]
      else: continue
    else:
      values = line.strip().split(args['separator'])
      sql = "%s " % action
      if action == "INSERT": sql += "INTO "
      sql += "`%s` SET " % table
      for i in range(len(fieldnames)):
        if fieldnames[i] == "NULL": continue
        if fieldnames[i] == "ID" and action == "UPDATE": continue
        if not sql.endswith(" "): sql += ", "
        sql += "%s='%s'" %(fieldnames[i], escape4sql(values[i]))
      if action == "UPDATE": sql += " WHERE ID='%s'" % values[fieldnames.index("ID")]
      print sql + ";"
  print "COMMIT;"

# =============================================================================
args = handle_arguments()
main( args )

