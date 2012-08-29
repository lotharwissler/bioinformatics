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
  stdout( " -f        nt alignment file" )
  stdout( " -t        tree file (newick format)" )
  stdout( " -p        path to PAML codeml" )
  stdout( " -n        number of cpus to use" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:t:p:n:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f':	args['aln'] = value
    if key == '-t':	args['tree'] = value
    if key == '-p':	args['codeml'] = value
    if key == '-n':	args['ncpu'] = int(value)
        
  if not args.has_key('aln'):
    stderr( "aln file missing." )
    show_help()
  if not file_exists( args.get('aln') ):
    stderr( "aln file does not exist." )
    show_help()
    
  if not args.has_key('tree'):
    stderr( "tree file missing." )
    show_help()
  if not file_exists( args.get('tree') ):
    stderr( "tree file does not exist." )
    show_help()

  if not args.has_key('codeml'):
    args['codeml'] =  "~/Results/Orthologs/PAML/codeml"
  args['pamlfolder'] = os.path.split(args.get('codeml'))[0] + '/'

  if not file_exists( args.get('codeml') ):
    stderr( "codeml binary not found." )
    show_help()
  if not dir_exists( args.get('pamlfolder') ):
    stderr( "paml folder does not exist" )
    show_help()

  return args

# =============================================================================
# =============================================================================
def main( args ):
  
  paml_orth_aln = args.get('pamlfolder')+'orth.aln'
  paml_orth_tree = args.get('pamlfolder')+'orth.tree'
  paml_orth_out = args.get('pamlfolder')+'orth.out'
  if file_exists(paml_orth_aln): os.unlink( paml_orth_aln )
  if file_exists(paml_orth_tree): os.unlink( paml_orth_tree )
  if file_exists(paml_orth_out): os.unlink( paml_orth_out )

  
  # copy all necessary files in the PAML folder
  #os.system( 'cp %s %s' % (args.get('aln'), paml_orth_aln) )
  os.system( 'sed s/^\>//g %s > %s' %( args.get('aln'), paml_orth_aln ) )
  os.system( 'cp %s %s' % (args.get('tree'), paml_orth_tree) )
  
  # now run PAML with a given model (ctl file) and return the results of the run
  models = ["M0", "Free", "M3K2", "M3K3", "M7", "M8"]
  #models = ["M3K2", "M3K3", "M7", "M8", "Free"]
  #models = ["M0"]
  #models = ["Free"]
  sys.stderr.write('%s\trunning PAML.codeml' % (os.path.split(args.get('aln'))[1]) )
  for M in models:
    sys.stderr.write('\t' + M)
    sys.stderr.flush()
    while not file_exists(paml_orth_out) or not file_exists(args.get('pamlfolder')+'rst'):
      CWD = os.getcwd()
      os.chdir( args.get('pamlfolder') )
      #os.system( 'cp %s codeml.ctl' % ('codeml.ctl.'+M) )
      error = os.WEXITSTATUS(os.system( './codeml codeml.ctl.' + M + ' &> codeml.log'))
      os.chdir( CWD )
    
    if not error:
      os.system( 'mv %s %s' % (paml_orth_out, args.get('aln')+'.paml.out.'+M) )
      os.system( 'mv %s %s' % (args.get('pamlfolder')+'rst', args.get('aln')+'.paml.rst.'+M) )
    else:
      os.system( 'mv %s %s' % (args.get('pamlfolder')+'codeml.log', args.get('aln')+'.paml.err.'+M) )
      sys.stderr.write(' (!)')
    sys.stderr.flush()
  sys.stderr.write("\n")

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
