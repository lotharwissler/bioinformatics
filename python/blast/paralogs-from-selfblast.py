#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
import tempfile
import hashlib
import fasta
from low import *			# custom functions, written by myself

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> -b <path>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        fasta file (input for blast)" )
	stdout( " -b        blast.out in tab format (-m 8)" )
	stdout( " -l        filter local (based on blast alignment)" )
	stdout( " -g        filter global (usign muscle)" )
	stdout( " -q        quiet: not stderr status messages" )
	stdout( " -i        add percent identity as a column" )
	stdout( " " )
	sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hf:b:qgli" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {'quiet':False, 'global':False, 'local':False, 'identity':False}
	for key, value in keys:
		if key == '-f': args['fasta'] = value
		if key == '-b': args['blastout'] = value
		if key == '-q': args['quiet'] = True
		if key == '-l': args['local'] = True
		if key == '-g': args['global'] = True
		if key == '-i': args['identity'] = True
				
	if not args.has_key('blastout'):
		stderr( "blastout file file missing." )
		show_help()
	if not file_exists( args.get('blastout') ):
		stderr( "blastout file file does not exist." )
		show_help()
		
	if not args.has_key('fasta'):
		stderr( "fasta file file missing." )
		show_help()
	if not file_exists( args.get('fasta') ):
		stderr( "fasta file file does not exist." )
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
class ParalogPair:
  def __init__(self, id1, id2):
    a = [id1, id2]
    a.sort()
    self.ids = a
    self.key = string.join(a, ",")
    self.pep = []
    self.name = hashlib.md5(string.join(self.ids, ",")).hexdigest() 
    
  def align_pep(self, overwrite=False):
    pf, pepfile = tempfile.mkstemp(".pep")
    af, alnfile = tempfile.mkstemp(".aln")
    os.close(pf)
    os.close(af)
    pf = open(pepfile, 'w')
    for i in range(len(self.pep)):
      pf.write(">" + self.ids[i] + "\n")
      pf.write(self.pep[i] + "\n")
    pf.close()
    os.system("muscle -in %s -out %s -quiet -maxiters 2 2> /dev/null" %(pepfile, alnfile))
    os.unlink(pepfile)
    self.aln = []
    af = open(alnfile)
    for line in af:
      if line.startswith(">"): self.aln.append("")
      else: self.aln[-1] += line.strip()
    af.close()
    os.unlink(alnfile)
    self.alnlen = len(self.aln[0])
    self.alnres = 0
    identity = 0
    for i in range(self.alnlen):
      x = self.aln[0][i]
      y = self.aln[1][i]
      if x != '-' and y != '-': 
        self.alnres += 1
        if x == y: identity += 1
    self.identity = 100.0 * identity / self.alnlen
    

# =============================================================================
def prefetch_sequences(pepfile):
  return fasta.get_sequence_hash(pepfile)

# =============================================================================
def get_seq_lengths(fastafile):
  lenhash = {}
  seqhash = fasta.get_sequence_hash(fastafile)
  for gid, seq in seqhash.iteritems(): lenhash[gid] = len(seq)
  return lenhash

# =============================================================================
def get_total_line_count(ifile):
  total = 0
  fo = open(args['blastout'])
  for line in fo: total += 1
  fo.close()
  return total

# =============================================================================
# =============================================================================
def main( args ):
  pephash, lenhash = prefetch_sequences(args['fasta']), get_seq_lengths(args['fasta'])
  current, total = 0, get_total_line_count(args['blastout'])
  fo = open(args['blastout'])
  for line in fo:
    current += 1
    if not args['quiet']: statusbar(current, total, "processing blastout")
    line = line.strip()
    if line.startswith("#") or len(line) == 0: continue
    (sid1, sid2, identity, alnlen, mismatch, gap, start1, stop1, start2, stop2, evalue, bitscore) = line.split("\t")
    if sid1 == sid2: continue
    if float(evalue) > float('1e-10'): continue
    if float(identity) < 30: continue
    
    pp = ParalogPair(sid1, sid2)
#    print >> sys.stderr, sid1 + "\t" + sid2
    if args['local']:
      length1, length2 = lenhash[sid1], lenhash[sid2]
      if length1 < 100 or length2 < 100: continue
      if int(alnlen) < 80 or (float(alnlen) / max([length1, length2])) < 0.70: continue
      #if alnlen < 100: continue
      pp.identity = float(identity)
      pp.alnres = float(bitscore)
    
    if args['global']:
      for pid in pp.ids: pp.pep.append(pephash[pid])
      pp.align_pep()
      if pp.alnres < 100 or pp.identity < 40: continue

    sys.stdout.write(string.join(pp.ids, "\t"))
    if args['identity']: print "\t" + str(pp.identity)
    else: print ""
    

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
