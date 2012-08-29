#!/usr/bin/python
# for all predicted peptides of a genome, compare its hmmout to the hmmout
# of all possible translations. putting this onto the scaffolds will
# show where genes might have been missed.

from low import *
import getopt, sys
import string
from gff3 import GeneFeature


# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hg:p:a:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-g': args['gff'] = value
    if key == '-p': args['predicted'] = value
    if key == '-a': args['alltranslations'] = value
    
  if not args.has_key('predicted'):
    stderr( "predicted hmmout file argument missing." )
    show_help()
  elif not file_exists( args.get('predicted') ):
    stderr( "predicted hmmout file does not exist." )
    show_help()

  if not args.has_key('alltranslations'):
    stderr( "alltranslations hmmout file argument missing." )
    show_help()
  elif not file_exists( args.get('alltranslations') ):
    stderr( "alltranslations hmmout file does not exist." )
    show_help()

  if not args.has_key('gff'):
    stderr( "gff file argument missing." )
    show_help()
  elif not file_exists( args.get('gff') ):
    stderr( "gff file does not exist." )
    show_help()


  return args

# =============================================================================
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -g <path> -p <path> -a <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -g        gff file of the genome under investigation" )
  stdout( " -p        hmmout of the predicted peptides" )
  stdout( " -a        hmmout of all translations" )
  stdout( " " )
  sys.exit(1)


# =============================================================================
def load_hmmout(file, gff=None):

  pid2domains = {}
  pid2scaffold = None
  if gff: pid2scaffold = load_gff(gff)
  fo = open(file, 'r')
  for line in fo:
    line = line.rstrip()
    columns = line.split()
    if len(columns) == 16: columns.pop(0) # remove first column (prot length)
    seqid, start, stop = columns[0], int(columns[3]), int(columns[4])

    if seqid.endswith("]"):
      posinfo = seqid[seqid.rindex("[")+1:-1]
      seqid = seqid[:seqid.rindex("[")]
      extractstart, extractstop = [int(e) for e in posinfo.split(":")[0:2]]
      if seqid == "scaffold_192": print "S:", extractstart, extractstop, start, stop
      if extractstop < extractstart: # indicating negative frame / anti-strand
        start = extractstop - (3*start)
        stop = extractstop - (3*stop)
        #stop, start = start, stop
      else:
        start = extractstart + (3*start)
        stop = extractstart + (3*stop)
      if seqid == "scaffold_192" and start > 10000: print "E:", extractstart, extractstop, start, stop

    elif pid2scaffold and pid2scaffold.has_key(seqid):
      exons = pid2scaffold[seqid]
      startnt = start * 3
      stopnt = stop * 3
      sumlength = 0
      for (scaffold, estart, estop) in exons:
        seqid = scaffold
        if not pid2domains.has_key(seqid): pid2domains[seqid] = []
        sumlength += (estop - estart)
        if startnt > sumlength: 
          #print scaffold, start, "(", startnt, ")", ">", sumlength
          continue
        if stopnt < sumlength - (estop - estart): 
          #print scaffold, stop, "(", stopnt, ")", "<", sumlength - (estop - estart)
          continue
        nstart = max([estop - (sumlength - startnt), estart])
        nstop  = min([estop, estop - (sumlength - stopnt)])
        pid2domains[seqid].append( [nstart, nstop] )

    if pid2scaffold: continue # already added the stuff
    if not pid2domains.has_key(seqid): pid2domains[seqid] = []
    pid2domains[seqid].append( [start, stop] )
  fo.close()
  return pid2domains

# =============================================================================
def load_gff(file):
  pid2scaffold = {}
  fo = open(file, 'r')
  for line in fo:
    if line.startswith("#"): continue
    feat = GeneFeature(line)
    #print >> sys.stderr, "feat type: %s" % feat.type
    if not feat.type == "exon": continue
    pid = feat.get_attributes()['Parent']
    pid = pid[pid.index(":")+1:]
    if not pid2scaffold.has_key(pid): pid2scaffold[pid] = []
    seqid, start, stop = feat.seqid, feat.start, feat.stop
    pid2scaffold[pid].append([seqid, start, stop])
  fo.close()
  print >> sys.stderr, "gff loaded with %s protein ids of relevance" % len(pid2scaffold)
  return pid2scaffold

# =============================================================================
def report_new(scaffold, hits):
  for (start, stop) in hits:
    print string.join([scaffold, str(start), str(stop)], "\t")

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  predictedHash = load_hmmout(args['predicted'], args['gff'])
  alltranslationsHash = load_hmmout(args['alltranslations'])
  print >> sys.stderr, "finished loading hmmouts."
  for scaffold, alltranshits in alltranslationsHash.iteritems():
    if not predictedHash.has_key(scaffold):
      report_new(scaffold, alltranshits)
      continue
    predhits = predictedHash[scaffold]
    print scaffold
    print alltranshits
    print ""
    print predhits
    sys.exit(99)


# =============================================================================
args = handle_arguments()
main( args )
