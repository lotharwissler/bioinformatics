#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself
import tempfile

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -a        nt alignment file" )
  stdout( " -s        evolver simulation file (generated with swapsc-in)" )
  stdout( " -t        tree file containing the newick tree (with numbers, not names, no spaces)" )
  stdout( " -p        path to swapsc binary" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "ha:s:p:t:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-a':	args['aln'] = value
    if key == '-s':	args['simul'] = value
    if key == '-p':	args['swapsc'] = value
    if key == '-t':	args['tree'] = value
        
  if not args.has_key('aln'):
    stderr( "aln file missing." )
    show_help()
  if not file_exists( args.get('aln') ):
    stderr( "aln file does not exist." )
    show_help()
    
  if not args.has_key('simul'):
    stderr( "simulation file missing." )
    show_help()
  if not file_exists( args.get('simul') ):
    stderr( "simulation file does not exist." )
    show_help()

  if not args.has_key('tree'):
    stderr( "tree file missing." )
    show_help()
  if not file_exists( args.get('tree') ):
    stderr( "tree file does not exist." )
    show_help()
  
  if not args.has_key('swapsc'):
    stderr("path to SWAPCS not specified")
    show_help()

  if not file_exists( args.get('swapsc') ):
    stderr( "SWAPSC binary not found." )
    show_help()

  args['workdir'] = os.path.split( args.get('swapsc') )[0] + '/'

  return args

# =============================================================================
def get_sysout(command):
  proc = os.popen(command)
  out = proc.read().strip()
  proc.close()
  return out 


# =============================================================================
def generate_control_file(outFolder):
  path = outFolder
  if not path.endswith("/"): path += "/"
  path += "SWAPSC.ctl"
  fw = open(path, "w")
  fw.write("""data_file: aln  *File with the alignment of sequences
  Tree_file: tree   *File with the phylogenetic tree in Newick format
  Output_file: out  *Name of the file with the output results
  Simulations: simul *Name of the file with the simulated alignments
  Model : 0       * 0 = Li 1993, 1 = Nei&Gojobori, 2 = Pamilo&Bianchi
  Window: 0     * 0 = inferred as in Fares et al. (2002), 1 = fixed
  Window_size: 3      * Size in codons of the window if fixed\n""")
  fw.close()


# =============================================================================
def main( args ):
  aln = args['aln']
  tree = args['tree']
  swapscbin = args['swapsc']
  simul = args['simul']
  outfile = aln + '.swapsc.out'

  tempdir = tempfile.mkdtemp(suffix='', prefix='tmp.swapsc.', dir='.')
  os.system("cp %s %s" %(aln, tempdir + '/aln'))
  os.system("cp %s %s" %(tree, tempdir + '/tree'))
  os.system("cp %s %s" %(simul, tempdir + '/simul'))
  generate_control_file(tempdir)
  os.chdir(tempdir)
  os.system(swapscbin + '&> swapsc.log')
  os.chdir("..")
  os.system("mv %s/out %s" %(tempdir, outfile))
  # check if output complete
  lastline = get_sysout("tail -n 1 %s" % outfile)
  if not lastline.startswith("P(neutral sites)"):
    os.system("mv %s/log %s.log" %(tempdir, outfile))
  os.system("rm -rf %s" % tempdir)


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
