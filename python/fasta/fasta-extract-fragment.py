#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file" )
  stdout( " -i        seq ID" )
  stdout( " -u        startpos (BLAST-like, starting at 1)" )
  stdout( " -v        endpos (BLAST-like)" )
  stdout( " -x        do not count gaps (BLAST-like)" )
  stdout( " -C        return complement (antisense strand sequence instead of sense)" )
  stdout( " -R        return reverse sequence" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """

  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:i:u:v:xCR" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'countGaps':True, 'complement':False, 'reverse':False}
  for key, value in keys:
    if key == '-f': args['fastafile'] = value
    if key == '-i': args['seqid'] = value
    if key == '-u': args['startpos'] = int(value) -1
    if key == '-v': args['endpos'] = int(value)
    if key == '-g': args['countGaps'] = False
    if key == '-C': args['complement'] = True
    if key == '-R': args['reverse'] = True
    
  for key in ['fastafile', 'seqid']:
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
def extract_sequence(args):
  fo = open(args['fastafile'])
  seqid = args['seqid']
  found = False
  for line in fo:
    if line.startswith('>'):
      line = line.rstrip()
      fid = line[1:].split()[0]
      if fid == seqid: found = True
      elif found: break
    if found: print line
  fo.close()
        
# =============================================================================
def extract_fragment(args):
  fo = open(args['fastafile'])
  seqid = args['seqid']
  pos = 0
  found = False
  startpos = args.get('startpos', 0)
  endpos = args.get('endpos', False)
  seq = ""
  for line in fo:
    line = line.rstrip()
    if line.startswith('>'):
      fid = line[1:].split()[0]
      if fid == seqid: 
        found = True
        out = '>' + fid
        if args.has_key('startpos') or args.has_key('endpos'): out += " %s:%s" %(startpos, endpos)
        if args['reverse']: out += " reverse"
        if args['complement']: out += " complement"
        print out
      elif found: break
      else: continue
    elif found: 
      if not args['countGaps']: line = line.replace('-','')
      if pos > endpos: break
      if pos < startpos and pos+len(line) < startpos: 
        pos += len(line)
        continue
      if pos < startpos and pos+len(line) >= startpos: out = line[startpos-pos:]
      elif pos >= startpos: out = line
      if not endpos == False and endpos < pos+len(line): out = out[:endpos-pos]
      if not args['reverse'] and not args['complement']:
        print out
      else:
        seq += out
      pos += len(line)
  fo.close()
  if args['reverse'] or args['complement']:
    seq = Seq(seq, IUPAC.unambiguous_dna)
    if args['reverse'] and args['complement']:
      seq = seq.reverse_complement()
    elif args['complement']: seq = seq.complement()
    elif args['reverse']: seq = seq[::-1]
    seq = str(seq)
    print seq

  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  if args.has_key('startpos') or args.has_key('endpos'):
    extract_fragment(args)
  else:
    extract_sequence(args)
  
  

# =============================================================================
args = handle_arguments()
main( args )

