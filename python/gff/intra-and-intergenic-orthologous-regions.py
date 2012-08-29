#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print >> sys.stderr, "usage: " + sys.argv[0] + " -f <gff-file> -o <orth-file>"
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        sorted all.parsed.gff file (species, chr, startpos)!!!" )
  stdout( " -o        clustered flybase gene orthologs file" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:o:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f': args['gff'] = value
    if key == '-o': args['orth'] = value
    
  if not args.has_key('gff'):
    print >> sys.stderr, "gff file argument missing."
    show_help()
  elif not file_exists( args.get('gff') ):
    print >> sys.stderr, "gff file does not exist."
    show_help()

  if not args.has_key('orth'):
    print >> sys.stderr, "orth file argument missing."
    show_help()
  elif not file_exists( args.get('orth') ):
    print >> sys.stderr, "orth file does not exist."
    show_help()

  return args

# =============================================================================
class Gene:
  def __init__(self, line):
    cols = line.rstrip().split("\t")
    self.species, self.name, self.chr, self.start, self.stop, self.strand = cols[0:6]
    self.start, self.stop = int(self.start), int(self.stop)
    self.prev, self.next = 0, 0
    self.orthologs = []

  def __cmp__(self, other):
    return cmp(self.start, other.start)

  def is_orthologous_to(self, other):
    if self in other.orthologs and other in self.orthologs: return 1
    return 0


# =============================================================================
def parse_gene_order(file):
  name2gene = {}
  fo = open(file)
  prevgene, prevspecies, prevchr = 0,0,0
  for line in fo:
    g = Gene(line)
    name2gene[g.name] = g
    if g.species == prevspecies and g.chr == prevchr:
      g.prev = prevgene
      prevgene.next = g
    prevgene, prevspecies, prevchr = g, g.species, g.chr
  fo.close()
  return name2gene

# =============================================================================
def get_orthologs(file, name2gene):
  fo = open(file)
  for line in fo: 
    if line.startswith("#"): continue
    if len(line.rstrip()) == 0: continue
    columns = line.rstrip().split("\t")
    genenames = [e[:e.index("(")] for e in columns]
    for gn in genenames:
      name2gene[gn].orthologs = [name2gene[x] for x in genenames]
      name2gene[gn].orthologs.remove(name2gene[gn])
  fo.close()
  return name2gene


# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):

  def process_gff_line(line, species):
    if line.startswith("#") or len(line.rstrip()) == 0: return
    columns = line.rstrip().split("\t")
    if len(columns) != 9: return
    type = columns[2]
    if type != "gene": return
    chr, start, stop, strand, descr = columns[0], columns[3], columns[4], columns[6], columns[8]
    id = re.search("ID=([^;]+);", descr).group(1)
    sys.stdout.write(species + "\t" + id + "\t")
    print string.join([chr, start, stop, strand], "\t")

# =============================================================================
    
  name2gene = parse_gene_order(args['gff'])
  name2gene = get_orthologs(args['orth'], name2gene)
  caught = {}
  for qname, qgene in name2gene.iteritems():
    # continue if already caught, or no neighbor
    if caught.has_key(qname): continue
    if not qgene.next: continue
    ngene = qgene.next
    # now check all direct orthologs of the query gene and see whether their neighbor is orthologous to the query's neighbor
    for ogene in qgene.orthologs:
      if not ogene.next: continue
      ongene = ogene.next
      if ngene.is_orthologous_to(ongene):
        print string.join([qgene.species, ogene.species, "intergenic", qgene.chr, str(qgene.stop +1), str(ngene.start -1), ogene.chr, str(ogene.stop +1), str(ongene.start -1)], "\t")
        print string.join([qgene.species, ogene.species, "gene", qgene.chr, str(qgene.start), str(qgene.stop), ogene.chr, str(ogene.start), str(ogene.stop)], "\t")
        print string.join([ngene.species, ongene.species, "gene", ngene.chr, str(ngene.start), str(ngene.stop), ongene.chr, str(ongene.start), str(ongene.stop)], "\t")
        caught[qgene.name] = 1
        caught[ngene.name] = 1
        caught[ogene.name] = 1
        caught[ongene.name] = 1
       
# =============================================================================
args = handle_arguments()
main( args )

