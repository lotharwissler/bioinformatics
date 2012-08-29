#!/usr/bin/python
# for an assembled genome (contigs/scaffolds/chromosomes), get all possible
# translations longer than a user-specified threshold.

from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate, Seq
from Bio.Alphabet import IUPAC
from low import *
import getopt, sys
import string


# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:t:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'minlength':30}
  for key, value in keys:
    if key == '-f': args['fasta'] = value
    if key == '-t': args['minlength'] = int(value)
    
  if not args.has_key('fasta'):
    stderr( "fasta file argument missing." )
    show_help()
  elif not file_exists( args.get('fasta') ):
    stderr( "fasta file does not exist." )
    show_help()

  return args

# =============================================================================
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> [-t <n>]" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        fasta file of the nucleotide sequences from which to predict the translations" )
  stdout( " -t        minimum length of a translation (amino acids) to be reported [default: 30]" )
  stdout( " " )
  sys.exit(1)


def process_seq(header, seq):
  hits = 0
  id = header
  if id.count(" ") > 0: id = id[:id.index(" ")]
  seq = Seq(seq)
  # direction1 is the direction we originally have had, 2 is the antisense strand
  # then TRANSLATE ALL POSSIBLE ORFs, do not stop at STOP codons
  dna_sequence_direction1 = seq
  dna_sequence_direction2 = dna_sequence_direction1.reverse_complement()
  translations = {}
  translations['+1'] = translate(dna_sequence_direction1)
  translations['-1'] = translate(dna_sequence_direction2)
  translations['+2'] = translate(dna_sequence_direction1[1:])
  translations['-2'] = translate(dna_sequence_direction2[1:])
  translations['+3'] = translate(dna_sequence_direction1[2:])
  translations['-3'] = translate(dna_sequence_direction2[2:])
  # get all polypeptides between stops, filter out those shorter than minlength
  polypeptides = {}
  for frame, translation in translations.iteritems():
    peptides = translation.split('*')
    if int(frame) < 0: startpos = len(seq) +1 + int(frame)
    else: startpos = int(frame)
    #print >> sys.stderr, "frame: %s | startpos: %s | scaffold length: %s" %(frame, startpos, len(seq))
    #print >> sys.stderr, "# peptides: %s | pep.length: %s | transformed length: %s | scaffold length: %s" %(len(peptides), sum([len(pep) for pep in peptides]), (sum([len(pep) for pep in peptides])+len(peptides))*3, len(seq))
    for peptide in peptides:
      peptide += '*'
      if int(frame) < 0: stoppos = startpos +1 - (3*len(peptide))
      else: stoppos = startpos -1 + (3*len(peptide))
      polypeptides[str(startpos)+':'+str(stoppos)] = peptide.tostring()
      if int(frame) < 0: startpos = stoppos-1
      else: startpos = stoppos+1

  for key, pepseq in polypeptides.iteritems():
    if len(pepseq) < args['minlength']: continue
    startpos, stoppos = [int(e) for e in key.split(":")]
    hits += 1
    print ">%s[%s:%s]" %( id, startpos, stoppos )
    print pepseq
  return hits

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  
  npep = 0
  fo = open(args['fasta'], 'r')
  id, seq = '', ''
  for line in fo:
    line = line.rstrip()
    if line.startswith(">"):
      if id != '' and seq != '': npep += process_seq(id, seq)
      id = line[1:]
      seq = ''
      sys.stderr.write("\r\tpeptides caught: %s " % npep)
    else: seq += line.strip()
  if id != '' and seq != '': npep += process_seq(id, seq)
  sys.stderr.write("\r\tpeptides caught: %s \n" % npep)

# =============================================================================
args = handle_arguments()
main( args )
