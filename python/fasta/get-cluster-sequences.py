#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself
import anydbm

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        orth.file" )
  stdout( " -d        dbm file (indexed fasta file)" )
  stdout( " -e        file extension of the fasta files to create (e.g. '.aa')" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:d:e:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f':	args['file'] = value
    if key == '-d':	args['dbm'] = value
    if key == '-e':	args['ext'] = value
        
  if not args.has_key('file'):
    stderr( "cluster file missing." )
    show_help()
  if not file_exists( args.get('file') ):
    stderr( "cluster file does not exist." )
    show_help()
    
  if not args.has_key('dbm'):
    stderr( "dbm file missing." )
    show_help()
  if not file_exists( args.get('dbm') ):
    stderr( "dbm file does not exist." )
    show_help()

  if not args.has_key('ext'):
    args['ext'] = '.fasta'


  return args

#
def get_shorthand( id ):
  if id.startswith("AT"):
    return "Arath"
  elif id.startswith("Pooc_"):
    return "Pooc"
  elif id.startswith("Zoma_"):
    return "Zoma"
  elif id.startswith("LOC_"):
    return "Oryza"
  elif id.startswith("IMGA|"):
    return "Medicago"
  elif id.startswith("jgi|Phypa1_1|"):
    return "Physco"
  elif id.startswith("jgi|Poptr1_1|"):
    return "Populus"
  elif id.startswith("GSVIV"):
    return "Vitis"
  elif id.startswith("AC"):
    return "Zea"
  elif id.startswith("Sbi"):
    return "Sorghum"
  else:
    stderr( "did not find a suitable shorthand for this id: %s" % id )

# =============================================================================
# =============================================================================
def main( args ):
  dbm = anydbm.open( args.get('dbm'), 'r' )
  fo = open( args.get('file') )
  counter = 0
  stdouttext, stderrtext = catch_bash_cmd_output("wc -l %s" %args.get('file') )
  #print "   stdout: %s" % stdouttext
  #print "   stderr: %s" % stderrtext
  clusternum = int(stdouttext.split()[0])

  for line in fo:
    counter += 1
    line = line.rstrip()
    ids = line.split("\t")[1:]
    c = add_leading_zeroes(counter, len(str(clusternum)))
    filename1 =  "orth.cluster.%s%s" % (c, args.get('ext'))
    filename2 =  "orth.cluster.%s%s" % (c, '.ids')
    fw = open( filename1, "w" )
    fwid = open( filename2, "w" )
    for id in ids:
      if not dbm.has_key( id ):
        stderr( "%s: id %s not available in the sequence hash" % (filename1, id) )
        continue
      fwid.write("%s\n" % id )
      fw.write(">%s\n%s\n" %( get_shorthand(id), dbm.get(id) ) )
    fw.close()
    fwid.close()
  fo.close()
  dbm.close()

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
