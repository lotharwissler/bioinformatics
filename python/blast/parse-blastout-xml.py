#!/usr/bin/python

import os, sys, getopt, string
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Alphabet import IUPAC

#==============================================================================
def show_help():
  print """%s parses BLASTX XML output to STDOUT
  
  Options:
  -f:\tBLASTX output in XML format
  -n:\tnumber of best hits to be parsed (default: 1)
  -e:\tmaximum e-value to accept hits (default: 1e-5)

	What this program does:
	It takes the best hit's start and endposition from BLAST, applies it to the sequence in your query (e.g. the CAP3-output),
	and translates to the left resp. right from the start resp. end of your CAP3-output, until a Start-orStopcodon appears.
  """ % sys.argv[0]

  sys.exit(1)


# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    sys.stderr.write( "no arguments provided.\n" )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:n:e:" )
  except getopt.GetoptError:
    sys.stderr.write( "invalid arguments provided.\n" )
    show_help()

  args = {}
  args['numhits'] = 1
  args['evalue'] = float('1e-5')
  for key, value in keys:
    if key == '-f': args['blastfile'] = value
    if key == '-n': args['numhits'] = int(value)
    if key == '-e': args['evalue'] = float(value)
    
  if not args.has_key('blastfile'):
    sys.stderr.write( "blastx XML file argument missing.\n" )
    show_help()
  elif not os.path.exists( args.get('blastfile') ) or not os.path.isfile( args.get('blastfile') ):
    sys.stderr.write( "blastx XML file does not exist.\n" )
    show_help()

  return args


#==============================================================================
def main(args):
  #print "Working..."
  header = ['query', 'hit', 'frame', 'query_startpos', 'query_endpos', 'subject_startpos', 'subject_endpos', 'evalue', 'score']
  print '#', string.join(header, "\t")
  XML = open( args.get('blastfile') )
  blast_records = NCBIXML.parse(XML)

  for i in blast_records:
  #  print i.query
    count = 0
    while count < args.get('numhits'):
      count += 1
      hit = i.alignments.pop(0)
      hsp = hit.hsps[0]
      if hsp.expect > args.get('evalue'): break
#      print i.query, hit.title.split()[0], hsp.frame[0], hsp.query_start, hsp.query_start -1+ len(hsp.query)*3, hsp.sbjct_start, hsp.sbjct_start -1+ len(hsp.sbjct), hsp.expect, hsp.score
      print string.join([i.query, hit.title.split()[0], 
        str(hsp.frame[0]), 
        str(hsp.query_start),
        str(hsp.query_start -1+ len(hsp.query.replace('-', ''))*3), 
        str(hsp.sbjct_start), 
        str(hsp.sbjct_start -1+ len(hsp.sbjct)), 
        str(hsp.expect),
        str(hsp.score)], "\t")


# =============================================================================
args = handle_arguments()
main( args )
