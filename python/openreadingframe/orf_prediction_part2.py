#!/usr/bin/python
# ORF prediction for the sequences without homologs

from Bio import SeqIO
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from Bio.Alphabet import IUPAC
import getopt, sys
import string


# getopt reads in my input file
input_file = ""

#print "All arguments: ", sys.argv
 
shortOptions = 'hf:t:'
longOptions = ['help', 'filename=', 'threshold=']
 
#==============================================================================
def usage():
    print """
    %s 
    -h\t--help\tdisplay help.
    -f\t--filename\tpath to FASTA input file in which to predict ORFs
    -t\t--threshold\tminimum length of predicted ORF (in amino acids) to be considered
    """ % sys.argv[0]


#==============================================================================
def get_parameters():

  if len(sys.argv) == 1:
    usage()
    sys.exit()

  opts = []
  args = []
  try:
      opts, args = getopt.getopt(sys.argv[1:], shortOptions, longOptions)
  except getopt.GetoptError:
      print "ERR: At least one option is not available!"
      usage()
      sys.exit()
   
  for o, a in opts:
    if o == "--help" or o == "-h":
      print "HELP"
      usage()
    elif o == "--filename" or o == "-f":
#        print "Filename:", a
      input_file = a
    elif o == "--threshold" or o == "-t":
      threshold = int(a)
   
  for a in args:
      print "Additional argument, no option: ", a    
  #end of getopt stuff! now it becomes even more exciting!!!

  return input_file, threshold


def get_input_sequences(input_file):
  ids2seqs = {}
  for seq_record in SeqIO.parse(open(input_file), "fasta"):
    ids2seqs[seq_record.id] = seq_record.seq
  return ids2seqs


############################################now the real programme starts########################################################

input_file, threshold = get_parameters()
ids2seqs = get_input_sequences(input_file)

#header = ['id', 'frame', 'startpos', 'endpos', 'cds', 'protein', 'evidence']
#print "#" + string.join(header, "\t")

# iterate input sequences
for key, dna_sequence_direction1 in ids2seqs.iteritems():
  # direction1 is the direction we originally have had, 2 is the antisense strand
  dna_sequence_direction2 = dna_sequence_direction1.reverse_complement()
  
  # TRANSLATE ALL POSSIBLE ORFs, do not stop at STOP codons
  translations = {}
  translations['1'] = translate(dna_sequence_direction1)
  translations['-1'] = translate(dna_sequence_direction2)
  translations['2'] = translate(dna_sequence_direction1[1:])
  translations['-2'] = translate(dna_sequence_direction2[1:])
  translations['3'] = translate(dna_sequence_direction1[2:])
  translations['-3'] = translate(dna_sequence_direction2[2:])

  polypeptides = {}
  for frame, translation in translations.iteritems():
    peptides = translation.split('*')
    startpos = 0
    for peptide in peptides:
      polypeptides[peptide.tostring()] = [frame, startpos]
      startpos += len(peptide)+1

  # get longest ORF with startpos and frame
  peptides = polypeptides.keys()
  peptides.sort(key=len)
  longestpeptide = peptides[-1]
  frame, startpos = polypeptides[longestpeptide]

  if len(longestpeptide) < threshold: continue

  start_nt = startpos *3
  stop_nt = start_nt + ((len(longestpeptide)+1)*3)
  if frame.startswith('-'):
    cds = dna_sequence_direction2.tostring()
  else:
    cds = dna_sequence_direction1.tostring()
  cds = cds[start_nt:stop_nt+1]

  if frame.startswith('-'):
    start_nt, stop_nt = stop_nt, start_nt
  
  outlist = [key, frame, str(start_nt), str(stop_nt), cds, longestpeptide, "2"]
  print string.join(outlist, "\t")
