#!/usr/bin/python

import os, sys, getopt, string, math
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Alphabet import IUPAC

OUTFILEPART2 = 'tmp.orf.part2.fasta'

#==============================================================================
def show_help():
  print """%s uses parsed BLASTX output to determine ORF, cds, and 
  putative protein sequence.
  
  Options:
  -f:\tFASTA file with the input nucleotide sequences  output in XML format
  -b:\tparsed BLASTX output with the best hit for each input sequence
  """ % sys.argv[0]

  sys.exit(1)


# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    sys.stderr.write( "no arguments provided.\n" )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:b::" )
  except getopt.GetoptError:
    sys.stderr.write( "invalid arguments provided.\n" )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['fastafile'] = value
    if key == '-b': args['blastfile'] = value
    
  if not args.has_key('blastfile'):
    sys.stderr.write( "parsed blastx best hit file argument missing.\n" )
    show_help()
  elif not os.path.exists( args.get('blastfile') ) or not os.path.isfile( args.get('blastfile') ):
    sys.stderr.write( "parsed blastx best hit file does not exist.\n" )
    show_help()

  if not args.has_key('fastafile'):
    sys.stderr.write( "input fasta file argument missing.\n" )
    show_help()
  elif not os.path.exists( args.get('fastafile') ) or not os.path.isfile( args.get('fastafile') ):
    sys.stderr.write( "input fasta file does not exist.\n" )
    show_help()

  return args


#==============================================================================
def get_blast_measures(blastfile):
  fo = open(blastfile)
  header = fo.readline()[1:].strip()
  # blaaa   blast -m8 format, where to get frame from?
  fields = ['query', 'hit', 'ident', 'aln_length', 'mismatches', 'gaps', 'query_startpos', 'query_endpos', 'subject_startpos', 'subject_endpos', 'evalue', 'score']
  hash = {}

  for line in fo:
    if line.startswith("#"): continue
    line = line.rstrip()
    elements = line.split("\t")
    if hash.has_key(elements[0]): continue
    hash[elements[0]] = {} 
    for i in range(len(fields)):
      hash[elements[0]][fields[i]] = elements[i]
    # define frame
    if int(hash[elements[0]]['query_startpos']) > int(hash[elements[0]]['query_endpos']):
      startpos, endpos = hash[elements[0]]['query_endpos'], hash[elements[0]]['query_startpos']
      hash[elements[0]]['query_startpos'], hash[elements[0]]['query_endpos'] = startpos, endpos
      startpos = int(startpos)
      frame = -1
    else:
      startpos = int(hash[elements[0]]['query_startpos'])
      frame = 1
    hash[elements[0]]['frame'] = str(frame * (startpos % 3))
  fo.close()

  return hash

#==============================================================================
def get_prototype_cds(prop, seq):
  if prop['frame'].startswith('-'):
    prototypeseq = Seq(seq, IUPAC.unambiguous_dna).reverse_complement().tostring()
#   print "#++", prop['query_endpos'], prop['query_startpos'], len(seq)
    orfstart = len(seq) - int(prop['query_endpos'])
    orfstop = len(seq) - int(prop['query_startpos'])
  else:
    prototypeseq = seq
    orfstart, orfstop = int(prop['query_startpos']) -1, int(prop['query_endpos'])

  return prototypeseq, orfstart, orfstop
 
#==============================================================================
def process_sequence(id, seq, blasthash):
  if not blasthash.has_key(id):
    fw = open(OUTFILEPART2, 'a')
    fw.write(">%s\n%s\n" %(id, seq))
    fw.close()
    return

  prop = blasthash[id]
  # get prototype cds
  # it is based on the blast hit but as long as possible
  # we walk backwards at the start position to then get the
  # sequence to analyze (sense if frame > 0, antisense if frame < 0)
  # and the adjusted start and stop positions of the minimum (homolog) orf
  prototypeseq, orfstart, orfstop = get_prototype_cds(prop, seq)
  #print "###", id, orfstart, orfstop
  tmpstart = orfstart
  while tmpstart > 3: tmpstart -= 3
  prototypeprotein = Seq(prototypeseq[tmpstart:], IUPAC.unambiguous_dna).translate().tostring()

  # determine start of ORF
  protstart = (orfstart - tmpstart) / 3

  #print >> sys.stderr, "####", 'id:', id, 'fr:', prop['frame'], 'qs:', prop['query_startpos'], 'qe:', prop['query_endpos'], 'O:', orfstart, 'P:', protstart, 'l:', len(prototypeprotein)# , prototypeprotein

  currentaminoacid = prototypeprotein[protstart]
  while protstart > 0 and currentaminoacid != 'M' and currentaminoacid != '*':
    protstart -= 1
    orfstart -= 3
    currentaminoacid = prototypeprotein[protstart]
    #print "#..", orfstart, currentaminoacid
  if currentaminoacid == '*': 
    protstart += 1
    orfstart += 3
    #print "#..", orfstart, prototypeprotein[protstart]

  # deterine end of ORF
  protstop = (orfstop - tmpstart) / 3
  if protstop >= len(prototypeprotein):
    protstop = len(prototypeprotein)-1
  #print "####", id, orfstart, orfstop, len(prototypeseq), protstop, len(prototypeprotein)
  currentaminoacid = prototypeprotein[protstop]
  while protstop < len(prototypeprotein) -1 and currentaminoacid != '*':
    protstop += 1
    orfstop += 3
    currentaminoacid = prototypeprotein[protstop]

  #cds = prototypeseq[orfstart:orfstop+1]
  if currentaminoacid == '*': protstop -= 1
  #protein = Seq(cds, IUPAC.ambiguous_dna).translate(to_stop=True).tostring()
  protein = prototypeprotein[protstart:protstop+1]
  cds = prototypeseq[orfstart:orfstart+(len(protein)*3)]
  
  #cds = prototypeseq[orfstart:]
  #seqtotranslate = prototypeseq[orfstart:]
  #protein = Seq(seqtotranslate, IUPAC.ambiguous_dna).translate(to_stop=False).tostring()
  #orfstop = orfstart + len(protein)*3
  #print "####", orfstart, orfstop, len(cds), len(protein)*3, protein
  #if protein.count('*') > 0:
  #  protein = protein[:protein.index('*')+1]
  #  orfstop = orfstart + len(protein)*3
  #  protein = protein[:-1]

  #print "####", orfstart, orfstop, len(seqtotranslate), len(protein)*3, protein

  #currentaminoacid = prototypeprotein[protstart]
  #cds = prototypeseq[orfstart:orfstop]
  # backtranslate positions for negative frames
  if prop['frame'].startswith('-'):
    startpos = len(prototypeseq) - orfstart
    endpos = len(prototypeseq) - orfstop
  else:
    startpos, endpos = orfstart, orfstop

  print string.join([id, prop['frame'], str(startpos), str(endpos), cds, protein, "1"], "\t")
  if len(protein) < 50:
    sys.stderr.write("protein shorter than 50aa!! check record %s (len=%s, %s)\n" %(id, len(protein), prop['frame'] ))

  if len(cds) < int(prop['aln_length']):
    print >> sys.stderr, "WARNING: cds shorter than NCBI blast hit. Check entry %s" % id
#  print "#####", len(cds), len(protein)*3, "\n"


#==============================================================================
def main(args):
  blasthash = get_blast_measures( args.get('blastfile') )
  #print blasthash
  header = ['record_name', 'frame', 'startpos', 'endpos', 'cds', 'protein', 'evidence']
  print '#', string.join(header, "\t")

  id, seq = "", ""
  fo = open( args.get('fastafile') )
  for line in fo:
    line = line.strip()
    if line.startswith('>'):
      if id != "": 
        process_sequence(id, seq, blasthash)
      
      id = line.split()[0][1:]
      seq = ""
    else: seq += line

  fo.close()
  if id != "": 
    process_sequence(id, seq, blasthash)

# =============================================================================
args = handle_arguments()
main( args )
