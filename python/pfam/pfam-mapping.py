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
  stdout( "usage: " + sys.argv[0] + " -f <path> " )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        pfam_full file to parse" )
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
    if key == '-f': args['file'] = value
    
  if not args.has_key('file'):
    stderr( "import file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "import file does not exist." )
    show_help()
  
  return args


# =============================================================================
def get_regex( args ):
  idhash = {}
  idhash['name'] = re.compile('^#=GF ID\s+(\S+)')
  idhash['acc'] = re.compile('^#=GF AC\s+(PF\S+)')
  idhash['descr'] = re.compile('^#=GF DE\s+(.*)$')
  idhash['comment'] = re.compile('^#=GF CC\s+(.*)$')
  idhash['pftype'] = re.compile('^#=GF TP\s+(\S+)')
  idhash['terminate'] = re.compile('^\\\\$')
  return idhash

# =============================================================================
class PfamEntry:
  def __init__(self):
    self.name = " "
    self.acc = " "
    self.descr = " "
    self.comment = " "
    self.pftype = " "

  def set_name(self, name):
    self.name = name

  def set_acc(self, acc):
    self.acc = acc

  def set_descr(self, descr):
    if self.descr == " ": self.descr = descr
    else: self.descr += " " + descr

  def set_comment(self, comment):
    if self.comment == " ": self.comment = comment
    else: self.comment += " " + comment

  def set_pftype(self, pftype):
    self.pftype = pftype
 
  def to_string(self):
    return string.join([self.name, self.acc, self.pftype, self.descr, self.comment], "\t")


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  regex = get_regex( args )
  info("getting total number of entries...")
  sout, serr = catch_bash_cmd_output( "grep '^//' -c %s" %args.get('file') )
  total = int( sout )
  info("total number of entries: %s\n" %total)
  count = 0
  infomsg("entries found: %s" %total )
  print "#" + string.join(['abbrv', 'pfamid', 'pfamtype', 'description', 'comment'], "\t")
  
  fo = open( args.get('file') )
  entry = PfamEntry()
  for line in fo:
    line = line.rstrip()
    if line.startswith('//'):
      count += 1
      print entry.to_string()
      info( "  status:\t%01.2f%%\t%s/%s" %( 100.0*count/total, count, total ) )
      entry = PfamEntry()
    elif not line.startswith('#=GF'): continue
    elif re.match(regex['name'], line): entry.set_name( re.match(regex['name'], line).group(1) )
    elif re.match(regex['acc'], line): entry.set_acc( re.match(regex['acc'], line).group(1) )
    elif re.match(regex['pftype'], line): entry.set_pftype( re.match(regex['pftype'], line).group(1) )
    elif re.match(regex['descr'], line): entry.set_descr( re.match(regex['descr'], line).group(1) )
    elif re.match(regex['comment'], line): entry.set_comment( re.match(regex['comment'], line).group(1) )
  fo.close()
  info("\ndone.\n")

# =============================================================================
args = handle_arguments()
main( args )
