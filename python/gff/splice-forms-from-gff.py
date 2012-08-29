#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print sys.stderr, "parses a gff file to map mRNA to genes."
  print sys.stderr, "depending on whether or not -p is given, it ouputs different results:"
  print sys.stderr, "-p absent:   reports the transcripts IDs for a given gene ID.  output: <geneid>tab<transcript1> <transcript2> <transcript3>..."
  print sys.stderr, "-p present:  reports the longest transcripts ID for a given gene ID.  output: <geneid>tab<transcript>\n"
  print sys.stderr, "usage: " + sys.argv[0] + " -f <gff-file> -p <peptides.fasta> [-n]"
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        gff file to parse" )
  stdout( " -p        peptide fasta file from which to extract the longest sequence for a gene with splice variants" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:p:n" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'name':0}
  for key, value in keys:
    if key == '-f': args['file'] = value
    if key == '-p': args['pep'] = value
    if key == '-n': args['name'] = 1
    
  if not args.has_key('file'):
    stderr( "gff file argument missing." )
    show_help()
  elif not file_exists( args.get('file') ):
    stderr( "gff file does not exist." )
    show_help()

  if args.has_key('pep') and not file_exists( args.get('pep') ):
    stderr( "peptide fasta file does not exist." )
    show_help()
    
  return args


# =============================================================================
def get_seq_lengths(file):
  lengthHash, id = {}, ""
  fo = open(file)
  for line in fo: 
    line = line.strip()
    if line.startswith(">"):
      id = line[1:]
      if id.count(" ") > 0: id = id[:id.index(" ")]
      if id.count("\t") > 0: id = id[:id.index("\t")]
      lengthHash[id] = 0 
    else: lengthHash[id] += len(line)
  return lengthHash



# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  def process_gff_line(line):
    columns = line.rstrip().split("\t")
    type, descrline = columns[2], columns[8]
    if type != "gene" and type != "mRNA" and type != "CDS": return
    currentdescr = {}
    for pair in descrline.split(";"): currentdescr[pair.split("=")[0]] = pair.split("=")[1]
    if type == "CDS":
      mRNAwithCDS[currentdescr["Parent"]] = 1
      return
    if not currentdescr.has_key("ID"):
      print >> sys.stderr, "ERROR: ID tag missing! line: \"" + line + "\""
      return
    if type == "mRNA" and not currentdescr.has_key("Parent"):
      print >> sys.stderr, "ERROR: Parent association missing! line: \"" + line + "\""
      return
    if type == "mRNA": gene2transcripts[currentdescr["Parent"]].append(currentdescr["ID"])
    #if currentdescr.has_key("Alias"): aliases[currentdescr["ID"]] = currentdescr["Alias"]

# =============================================================================

  aliases = {}
  gene2transcripts = defaultdict(list)
  mRNAwithCDS = {}
  fo = open( args.get('file') )
  for line in fo:
    if line.startswith("#"): continue
    process_gff_line(line)
  fo.close()

  if args.has_key('pep'):
    lengthHash = get_seq_lengths(args['pep'])
    for gene, associds in gene2transcripts.iteritems():
      galias = gene
      if aliases.has_key(gene): galias = aliases[gene]
      if len(associds) == 1:
        talias = associds[0]
        if aliases.has_key(talias): talias = aliases[talias]
        if not mRNAwithCDS.has_key(associds[0]): continue
        if not lengthHash.has_key(talias): talias = talias.split(":")[1]
        if not lengthHash.has_key(talias): talias = talias.split("-")[0]
        if not lengthHash.has_key(talias): 
          print >> sys.stderr, "ERROR: could not find fasta sequence in peptide file with ID \"" + talias + "\""
          continue
        print galias + "\t" + talias
        continue
      peptides = {}
      for associd in associds:
        talias = associd
        if aliases.has_key(talias): talias = aliases[talias]
        if not mRNAwithCDS.has_key(associd): continue
        if not lengthHash.has_key(talias): talias = talias.split(":")[1]
        if not lengthHash.has_key(talias): talias = talias.split("-")[0]
        if not lengthHash.has_key(talias): 
          print >> sys.stderr, "ERROR: could not find fasta sequence in peptide file with ID \"" + talias + "\""
          continue
        peptides[talias] = lengthHash[talias]
      if len(peptides) > 0:
        best = sorted(peptides.iteritems(), key=lambda x: x[1] , reverse=True)[0][0]
        print galias + "\t" + best
    


  else:
    for gene, associds in gene2transcripts.iteritems():
      galias = gene
      if aliases.has_key(gene): galias = aliases[gene]
      sys.stdout.write(galias + "\t")
      for i in range(len(associds)):
        talias = associds[i]
        if aliases.has_key(talias): 
          associds[i] = aliases[talias]
      sys.stdout.write(string.join(associds, " ") + "\n")




  

# =============================================================================
args = handle_arguments()
main( args )

