#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import sqlite3
from low import *  # custom functions, written by myself
from sciroko import SSR
from collections import defaultdict


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        sciroko output file" )
  stdout( " -d        db file" )
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
    if key == '-f': args['file'] = value
    if key == '-d': args['db'] = value
    
  if not args.has_key('file'):
    stderr( "sciroko file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "sciroko file does not exist." )
    show_help()
  
  return args

# =============================================================================
def init_db(conn):
  conn.execute("CREATE TABLE IF NOT EXISTS ssrs(id INTEGER PRIMARY KEY ASC, organism VARCHAR(50), chr VARCHAR(50), startpos INTEGER, endpos INTEGER, motif VARCHAR(10), motif_std VARCHAR(10), length INTEGER, score INTEGER, mismatches INTEGER, binpos INTEGER, seq VARCHAR(255))")


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  conn = sqlite3.connect(args['db'])
  init_db(conn)
  
  fo = open(args['file'])
  for line in fo:
    if not line.startswith("RAL"): continue
    m = SSR(line)
    sql = "INSERT INTO ssrs(organism, chr, startpos, endpos, motif, motif_std, length, score, mismatches, binpos, seq) VALUES (\'%s\', \'%s\', %s, %s, \'%s\', \'%s\', %s, %s, %s, %s, \'%s\')" %(m.organism, m.chromosome, m.startpos, m.endpos, m.motif, m.motif_std, m.length, m.score, m.mismatches, m.megabase, m.seq)
    conn.execute(sql)
  res = conn.execute("SELECT COUNT(*) FROM ssrs")
  entries = res.fetchall()[0][0]
  print "done. entries added:", entries
  conn.commit()
  conn.close()


# =============================================================================
args = handle_arguments()
main( args )

