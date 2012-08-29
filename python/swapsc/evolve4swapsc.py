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
  stdout( " -m        paml M0 out file" )
  stdout( " -p        path to PAML evolver" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:m:p:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f':	args['aln'] = value
    if key == '-m':	args['m0'] = value
    if key == '-p':	args['evolver'] = value
        
  if not args.has_key('aln'):
    stderr( "aln file missing." )
    show_help()
  if not file_exists( args.get('aln') ):
    stderr( "aln file does not exist." )
    show_help()
    
  if not args.has_key('m0'):
    stderr( "M0 file missing." )
    show_help()
  if not file_exists( args.get('m0') ):
    stderr( "M0 file does not exist." )
    show_help()

  if not args.has_key('evolver'):
    stderr("path to PAML evolver not specified")
    show_help()

  if not file_exists( args.get('evolver') ):
    stderr( "evolver binary not found." )
    show_help()

  return args

# =============================================================================
def generate_evolver_infile( file, settings ):
  setnames = ['mode', 'random seed number', '# sequences', '# sites', '# replicates', 'tree length', 'newick tree', 'omega', 'kappa', 'codons']
  for name in setnames:
    if not settings.has_key(name):
      stderr( 'key %s missing in settings hash' %(name) )
  fw = open( file, 'w' )
  line = [settings.get('mode'), "* mode", "\n"]
  fw.write(string.join(line, ' '))
  line = [settings.get('random seed number'), "* random seed number", "\n"]
  fw.write(string.join(line, ' '))
  line = [settings.get('# sequences'), settings.get('# sites'), settings.get('# replicates'), "* #seq  #sites  #replicates", "\n"]
  fw.write(string.join(line, ' '))
  fw.write("\n")
  line = [settings.get('tree length'), "* tree length", "\n"]
  fw.write(string.join(line, ' '))
  line = [settings.get('newick tree'), "\n"]
  fw.write(string.join(line, ' '))
  fw.write("\n")
  line = [settings.get('omega'), "* omega", "\n"]
  fw.write(string.join(line, ' '))
  line = [settings.get('kappa'), "* kappa", "\n"]
  fw.write(string.join(line, ' '))
  fw.write("\n")
  fw.write(settings.get('codons'))
  fw.write("\n")
  fw.write('// end of file.')
  fw.flush()
  fw.close()

# =============================================================================
def generate_tree_file( file, settings ):
  tree = settings.get('tree for file')
#  pos = tree.index('),(')+1
#  if tree.index('1') > pos:
#    part1 = tree[1:pos]
#    part2 = tree[pos+1:-1]
#    tree = '(' + part2 + ',' + part1 + ')' + "\n"

#  while re.search('\),\d+\)', tree):
#    print "rearranging tree: %s" % tree
#    match = re.search('\),(\d+)\)', tree)
#    xpos = match.start(1) -1
#    node = match.group(1)
#    print "pos = %s, node = %s" %(xpos,node)
#    pos = xpos
#    count = 0
#    while 1:
#      pos -= 1
#      if tree[pos] == ')': 
#        count += 1
#        continue
#      if tree[pos] == '(':
#        count -= 1
#      if count == 0:
#        tree = tree[:pos] + node + ',' + tree[pos:xpos-1] + ')' + tree[xpos+2:] + "\n"
#        break
  
  tree = tree.rstrip()
  #beginbrackets = len(re.search("^(\(+)", tree).group(1))
  #endbrackets = len(re.search("(\)+)$", tree).group(1))
  #if beginbrackets > endbrackets:
  #  xpos = tree.rindex("),(") +1
  #  tree = '(' + tree[xpos+1:-1] + ',' + tree[1:xpos] + ')'

  fw = open( file, "w" )
  fw.write( tree )
  fw.close()
  #print "final tree: %s" % tree

# =============================================================================
# =============================================================================
def main( args ):

  #sys.stderr.write(args.get('aln') + "\t")
  #sys.stderr.flush()
  # create evolver control file based on the M0 out file
  TARGET = args.get('aln')+'.evolver.out'
  if os.path.exists(TARGET) and os.path.isfile(TARGET) and os.path.getsize(TARGET) > 0: return

  settings = {'mode':'0', '# replicates':'1000'}
  fo = open( args.get('m0') )
  line = ""
  while not line.startswith('Time used:'):
    line = fo.readline()
    if line.startswith('seed used ='):
      settings['random seed number'] = re.match('seed used =\s*(\d+)', line).group(1)
      #line = fo.readline()
      while not re.match("\s+\d+\s+\d+\s*$", line): 
        line = fo.readline()
      numbers = line.split()
      settings['# sequences'], settings['# sites'] = numbers[0:2] 
      settings['# sites'] = str( (int(settings['# sites'])/3) )
      continue
    if line.startswith('Codon frequencies under model, for use in evolver'):
      line = fo.readline()
      settings['codons'] = ''
      while re.search("\S+", line):
        settings['codons'] += line
        line = fo.readline()
      continue
    if line.startswith('tree length ='):
      settings['tree length'] = re.match('tree length =\s+(\S+)', line).group(1)
      line = fo.readline()
      line = fo.readline()
      settings['newick tree'] = re.match('(.*)$', line).group(1)
      continue
    if line.startswith('TREE #  1:'):
      settings['tree for file'] = re.search('(\(.*\));', line).group(1).replace(' ','')
    if line.startswith('kappa'):
      settings['kappa'] = re.match('kappa\s+\(ts/tv\)\s+=\s+(\S+)', line).group(1)
      continue
    if line.startswith('omega'):
      settings['omega'] = re.match('omega\s+\(dN/dS\)\s+=\s+(\S+)', line).group(1)
      break
  fo.close()

  tmpFolder = "." + args.get('aln')
  if not os.path.exists( tmpFolder ): os.mkdir( tmpFolder )
  os.chdir( tmpFolder )
  
  evolver_in = 'evolver.in'
  generate_evolver_infile(evolver_in, settings)
  treefile = args.get('aln') + '.swapsc.tree'
  generate_tree_file(treefile, settings)
  #sys.exit(1)

  #sys.stderr.write('> evolver.in'+ "\t")
  #sys.stderr.flush()

  # run evolver, save output
  os.system( '%s 6 %s &> evolver.log' % (args.get('evolver'), evolver_in) )
  #sys.stderr.write('> evolver.out'+ "\t")
  #sys.stderr.flush()
  evolver_out = "mc.paml"
  if file_exists('evolver.out'): os.unlink('evolver.out')
  if file_exists('ancestral.txt'): os.unlink('ancestral.txt')
  if file_exists('evolver.log'): os.unlink('evolver.log')
  if file_exists('evolver.in'): os.unlink('evolver.in')
  if file_exists(evolver_out):
    os.system('mv %s %s' %(evolver_out, "../"+TARGET) )
    #sys.stderr.write('> success'+ "\n")
  #else:
    #sys.stderr.write('> ERROR'+ "\n")
  os.chdir("..")
  os.system(" rm -rf " + tmpFolder )


  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
