#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import * 

OUTPUTHASH = {
  'q'   : 'query',
  'h'   : 'hitid',
  'd'   : 'hitdescr',
  'sl'   : 'sbjct_length',
  'ql'   : 'query_length',
  'e'   : 'evalue',
  's'   : 'score',
  'qs'  : 'query_startpos',
  'qe'  : 'query_endpos',
  'ss'  : 'sbjct_startpos',
  'se'  : 'sbjct_endpos',
  'hl'   : 'hitlength',
  'i'   : 'identities',
  'p'   : 'positives',
  'g'   : 'gaps',
  'frm'  : 'frame',
  'str'  : 'strand',  
}


# =============================================================================  
def show_help():
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> [-e -n -i -p -l -o <string>] > out.file" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        blast.out file to be parsed" )
  stdout( " -n        number of best hits to parse" )
  stdout( " -e        minimum evalue of a hit to be parsed" )
  stdout( " -i        minimum identity (in %)" )
  stdout( " -l        minimum length of a hit to be parsed" )
  stdout( " -d        delimiter that is used in the stdout to seperate the fields. default: space" )
  stdout( "           also allowed: ';' and 'tab' and ',' and '|'" )
  stdout( " " )
  stdout( " -o        output fields. default: \"q,h,s,e,qs,qe,ss,se,hl,sl,i,p\"" )
  for k, v in OUTPUTHASH.iteritems(): stdout( "                          %s:\t%s" %(k,v) )
  stdout( " " )
  
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hb:f:n:e:l:i:p:o:d:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()
  
  args = {}
  for key, value in keys:
    if key == '-f':
      if not file_exists( value ):
        stderr( "invalid path in " + key )
        show_help()
      else:
        args['file'] = value
    if key == '-e':  args['evalthresh'] = float(value)
    if key == '-l':  args['minlength'] = int(value)
    if key == '-i':  args['minident'] = int(value)
    if key == '-p':  args['minpos'] = int(value)
    if key == '-n':  args['numberofbesthits'] = int(value)
    if key == '-o':  args['outformat'] = value
    if key == '-d':  args['delimiter'] = value
        
  if not args.has_key('file'):
    stderr( "blast out file missing." )
    show_help()
  
  if not args.has_key('outformat'):
    args['outformat'] = "q,h,s,e,qs,qe,ss,se,hl,sl,i,p"
  
  args['tooutput'] = []
  for key in args.get('outformat').split(','):
    if OUTPUTHASH.has_key(key):
      args['tooutput'].append( OUTPUTHASH.get(key) )
    else:
      stderr( 'OUTPUTHASH does not contain the key "%s"' %(key) )
  
  if not args.has_key('delimiter'):
    args['delimiter'] = " "
  
  return args

# =============================================================================
def print_hit(args,hithash):
  """
  returns counted (0,1)
  """
  
  if args.has_key('evalthresh'):
    if args.get('evalthresh') < float(hithash.get('evalue')):  return 0
  if args.has_key('minlength'):
    if args.get('minlength') > int(hithash.get('hitlength')): return 0
  if args.has_key('minident'):
    if args.get('minident') > int(hithash.get('identities')): return 0
  if args.has_key('minpos'):
    if args.get('minpos') > int(hithash.get('positives')): return 0
  
  L = []  
  for key in args.get('tooutput'):
    if hithash.has_key( key ):
      L.append(hithash.get( key ))
    else:
      stderr( 'hithash does not contain the key "%s"' %(key) )
      sys.exit(1)
  
  if args.get('delimiter') in [ 't', 'tab' ]:
    print string.join(L, '\t')
  elif args.get('delimiter') in [ ';' ',' '|' ]:
    print string.join(L, args.get('delimiter'))
  else:
    print string.join(L, ' ')
  
  return 1


# =============================================================================
class BlastOutput:
  def __init__(self, file):
    self.filehandle = open(file, 'r')

  def next_query(self):
    print "# next_query"
    self.querylines = []
    while 1:
      line = self.filehandle.readline()
      # find start
      if len(self.querylines) == 0 and not line.startswith("Query="): continue
      # break at end
      if line.startswith("BLAST") or line.startswith("  Database"): break
      # else append
      self.querylines.append(line.rstrip())

  def next_hit(self):
    print "# next_hit"
    self.hitlines = []
    index = 0
    while 1:
      if index == len(self.querylines): break
      line = self.querylines[index]
      index += 1
      # find start
      if len(self.hitlines) == 0 and not line.startswith(">"): continue
      # break at end
      if len(self.hitlines) > 0 and line.startswith(">"): break
      # else append
      self.hitlines.append(line.rstrip())
      print line.rstrip()

  def next_hsp(self):
    print "# next_hsp"
    self.hsplines = []
    index = 0
    while 1:
      if index == len(self.hitlines): break
      line = self.hitlines[index]
      index += 1
      # find start
      if len(self.hsplines) == 0 and not line.startswith(" Score ="): continue
      # break at end
      if len(self.hsplines) > 0 and line.startswith(" Score ="): break
      # else append
      self.hsplines.append(line.rstrip())


  def parse_querylines(self):
    queryid = self.querylines[0].split()[1]
    return queryid

  def parse_hitlines(self):
    elements = string.join(self.hitlines, " ")[1:].split()
    subject_id = elements[0]
    subject_length = elements[-1]
    subject_description = elements[1:-3]
    return subject_id, subject_description, subject_length

  def parse_hsplines(self):
    topdefs = string.join(self.hsplines[0:3], " ")
    hsp_frame = None
    hsp_strand = None
    hsp_score = re.search('Score =\S+(\d+)', topdefs).group(1)
    hsp_evalue = re.search('Expect[\(\)0-9]* =\S+([e\.0-9-]+)', topdefs).group(1)
    hsp_identity = re.search('Identities =\S+\d+/\d+\S+\((\d+)%\)', topdefs).group(1)
    hsp_positive = re.search('Positives =\S+\d+/\d+\S+\((\d+)%\)', topdefs).group(1)
    if re.search('Frame =', topdefs): hsp_frame = re.search('Frame =\S+([-+1-3]+)', topdefs).group(1)
    return hsp_score, hsp_evalue, hsp_identity, hsp_positive, hsp_frame, hsp_strand


  def parse(self):
    while 1:
      self.next_query()
      if not self.querylines[0].startswith("Query="): break
      while 1:
        self.next_hit()
        if not self.hitlines[0].startswith(">"): break
        while 1:
          self.next_hsp()
          if not self.hsplines[0].startswith(" Score ="): break
          queryid = self.parse_querylines()
          subject_id, subject_description, subject_length = self.parse_hitlines()
          hsp_score, hsp_evalue, hsp_identity, hsp_positive, hsp_frame, hsp_strand = self.parse_hsplines()
          print string.join([queryid, subject_id, hsp_score, hsp_evalue, hsp_frame], "\t")
      

    self.filehandle.close()


# =============================================================================
def parse_blast_out( args ):
  #print "# blast.out file:", args.get('file')
  #print "# numberofbesthits:", args.get('numberofbesthits')
  #print "# max.evalue:", args.get('evalthresh')
  #print "# min.length:", args.get('minlength')
  #print "# fields: query, hitid, score, evalue, query_startpos, query_endpos, sbjct_startpos, sbjct_endpos, hitlength, subjct_length, identities, positives, frame_or_strand"
  
  parser = BlastOutput( args.get('file') )
  parser.parse()
  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
parse_blast_out( args )
