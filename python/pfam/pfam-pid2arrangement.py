#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
import pfam
import copy


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        hmmout file" )
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
    if key == '-f': args['hmmoutfile'] = value
    
  for key in ['hmmoutfile']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key):
      print >> sys.stderr, "missing argument", key
      show_help()
  return args

# =============================================================================
def statusbar(current, total, message="", width=40):
  progress = 1.0*current/total
  if message != "": message = "[" + message + "]"
  progressbar = "=" * int(progress*width)
  while len(progressbar) < width: progressbar += " " 
  sys.stderr.write("\r   0% " + progressbar + " 100% " + message)
  if progress == 1.0: sys.stderr.write("\n")

# =============================================================================
def get_dom_coverage(domains):
  endpos = max([int(d.get_attr('alignment_end')) for d in domains])
  cov = [0] * endpos
  for d in domains:
    for i in range(int(d.get_attr('alignment_start')), int(d.get_attr('alignment_end'))):
      cov[i] += 1
  return cov
  
# =============================================================================
def domains2arrangement(domains):
  if len(domains) == 1: return domains[0].get_attr('hmm_name')
  domains.sort(cmp=lambda x,y: cmp(int(x.get_attr('alignment_start')), int(y.get_attr('alignment_start'))))
  domainCoverage = get_dom_coverage(domains)
  while max(copy.copy(domainCoverage)) > 1:
    pos = domainCoverage.index(max(copy.copy(domainCoverage)))
    resolveDomains = []
    for d in domains:
      if d.covers(pos): resolveDomains.append(d)
    resolveDomains.sort(cmp=lambda x,y: cmp(float(x.get_attr('E-value')), float(y.get_attr('E-value'))))
    domains.remove(resolveDomains[-1])
    domainCoverage = get_dom_coverage(domains)
  domarray = [d.get_attr('hmm_name') for d in domains]
  if len(domarray) > 1:
    i = 0
    while i < len(domarray)-1:
      if domarray[i] == domarray[i+1]: domarray.pop(i+1)
      else: i += 1
  return string.join(domarray, ';')    
  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  pid2domains = pfam.read_hmmout(args['hmmoutfile'])
  for pid, domains in pid2domains.iteritems():
    sys.stderr.write("\r" + pid + " "*20)
    arr = domains2arrangement(domains)
    print pid + "\t" + arr          

# =============================================================================
args = handle_arguments()
main( args )

