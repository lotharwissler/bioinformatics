#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import math        # match functions
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> -t -a -n -m [-i -p]" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file to import" )
  stdout( " -t        type: fasta, or xdom" )
  stdout( " -a        action: \"insert\" or \"update\"" )
  stdout( " -n        mysql table name" )
  stdout( " -m        field names (comma separated), mapping to the fields to be parsed in the same order, leave out ID" )
  stdout( " -p        prefix to put in fron of the key" )
  stdout( " -i        record name to id file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:n:t:m:i:p:a:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-a': args['action'] = value
    if key == '-t': args['type'] = value
    if key == '-n': args['table'] = value
    if key == '-m': args['fields'] = value
    if key == '-i': args['idfile'] = value
    if key == '-p': args['prefix'] = value
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
    
  if args.has_key('idfile') and not file_exists( args.get('idfile') ):
    stderr( "idfile does not exist." )
    show_help()
    

  if not args.has_key('type'):
    stderr( "file format argument missing." )
    show_help()
  
  if not args.has_key('action'):
    stderr( "action argument missing." )
    show_help()

  if not args.has_key('table'):
    stderr( "table name missing." )
    show_help()
   
  if not args.has_key('fields'):
    stderr( "field names missing." )
    show_help()

  return args


# =============================================================================
def get_ids( args ):
  """ 
  reads in the idfile and returns a hash, mapping record names to IDs
  """
  idhash = {}
  fo = open( args.get('idfile') )
  for line in fo:
    recordname, recordid = line.split()
    idhash[ recordname ] = recordid
  return idhash
  

# =============================================================================
def sql_out(action, table, fieldvaluelist):
  """ 
  writes an sql insert statement to stdout
  """
  count = 0

  if action == "insert":
    sys.stdout.write( "INSERT INTO `" + table + "` SET " )
    for kv in fieldvaluelist:
      if count > 0: sys.stdout.write(", ")
      sys.stdout.write( kv[0] + "='" + kv[1] + "'"  )
      count += 1

  elif action == "update":
    sys.stdout.write( "UPDATE `" + table + "` SET " )
    for kv in fieldvaluelist[1:]:
      if count > 0: sys.stdout.write(", ")
      sys.stdout.write( kv[0] + "='" + kv[1] + "'"  )
      count += 1
    sys.stdout.write( " WHERE " + fieldvaluelist[0][0] + "='" + fieldvaluelist[0][1] + "'" )

  sys.stdout.write(";\n")

# =============================================================================
def parse_fasta( args ):

  action = args.get('action')
  tablename = args.get('table')
  fields = args.get('fields').split(',')
  if not len(fields) == 2:
    stderr("expected 2 fields to parse fasta.")
    sys.exit(3)
  
  fo = open( args.get('file'), 'r' )
  key, value = "", ""
  for line in fo:
    line = line.rstrip().rstrip("\n")
    if line.startswith( ">" ):
      if key != "":
        sql_out(action, tablename, [ [fields[0], key], [fields[1], value] ])
        key, value = "", ""
      if args.has_key('prefix'): key = args.get('prefix') + line[1:]
      else: key = line[1:]
    else:
      value += line.strip()

  fo.close()

# =============================================================================
def parse_xdom( args ):
  """
  """
  action = args.get('action')
  tablename = args.get('table')
  fields = args.get('fields').split(',')
  if not len(fields) >= 2:
    stderr("expected 2 or more fields to parse fasta.")
    sys.exit(3)
  
  fo = open( args.get('file'), 'r' )
  fieldvaluelist = []
  for line in fo:
    line = line.rstrip().rstrip("\n")

    if line.startswith( ">" ):
      if args.has_key('prefix'): fieldvaluelist.append( [ fields[0], args.get('prefix') + line[1:] ] )
      else: fieldvaluelist.append( [ fields[0], line[1:] ] )

    else:
      values = line.split("\t")
      for i in range( len(fields[1:]) ):
        fieldvaluelist.append( [ fields[i+1], values[i] ] )
      sql_out(action, tablename, fieldvaluelist)
      fieldvaluelist = fieldvaluelist[:1]

  fo.close()
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  if args.has_key('idfile'): idhash = get_ids( args )
  fileformat = args.get('type')
  if fileformat.lower() == 'fasta':
    parse_fasta(args)
  elif fileformat.lower() == 'xdom':
    parse_xdom(args)
  else:
    stderr( "invalid file format. only fasta or xdom allowed." )
    sys.exit(2)


# =============================================================================
args = handle_arguments()
main( args )

