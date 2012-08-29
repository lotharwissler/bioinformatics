#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
import tempfile
from low import *			# custom functions, written by myself

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        nt alignment file" )
  stdout( " -t        tree file (newick format)" )
  stdout( " -m        models to run (comma separates)" )
  stdout( " -p        path to PAML codeml" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:t:m:p:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f':	args['aln'] = value
    if key == '-t':	args['tree'] = value
    if key == '-p':	args['codeml'] = value
    if key == '-m':	args['models'] = value.split(",")
        
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
  if not args.has_key('models'): 
    stderr( "no models to run." )
    show_help()

  if not file_exists( args.get('codeml') ):
    stderr( "codeml binary not found." )
    show_help()
  args['pamlfolder'] = os.path.split(args.get('codeml'))[0] + '/'
  if not dir_exists( args.get('pamlfolder') ):
    stderr( "paml folder does not exist" )
    show_help()

  return args

# =============================================================================
# =============================================================================
def main( args ):
  
  models = args['models']
  aln, tree = args['aln'], args['tree']
  codemlbin = args['codeml']
  ctlbase = args['pamlfolder'] + 'codeml.ctl.'
  for model in models:
    if tree.count("."):
      ext = os.path.splitext(tree)[1]
      outfile = aln+".codeml"+ext+"."+model
    else:
      outfile = aln+".codeml."+model
    tempdir = tempfile.mkdtemp(suffix=model, prefix='tmp.codeml.', dir='.')
    os.system("cp %s %s" %(aln, tempdir + '/in-aln'))
    os.system("cp %s %s" %(tree, tempdir + '/in-tree'))
    os.system("cp %s %s" %(ctlbase + model, tempdir + '/codeml.ctl'))
    os.chdir(tempdir)
    os.system(codemlbin)
    os.chdir("..")
    os.system("mv %s/out-codeml %s" %(tempdir, outfile))
    os.system("rm -rf %s" % tempdir)


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
