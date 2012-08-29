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
  for k, v in OUTPUTHASH.iteritems():
    stdout( "                          %s:\t%s" %(k,v) )
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
def parse_blast_out( args ):
  #print "# blast.out file:", args.get('file')
  #print "# numberofbesthits:", args.get('numberofbesthits')
  #print "# max.evalue:", args.get('evalthresh')
  #print "# min.length:", args.get('minlength')
  #print "# fields: query, hitid, score, evalue, query_startpos, query_endpos, sbjct_startpos, sbjct_endpos, hitlength, subjct_length, identities, positives, frame_or_strand"
  
  sout, serr = catch_bash_cmd_output( "grep 'Query=' -c %s" %args.get('file') )
  total = int( sout )
  count = 0
  hitcount = 0
  hithash = {}
  fh = open( args.get('file') )
  for line in fh:
    # new hit entry
    if (line.startswith('Query=') or line.startswith('>')) and len(hithash) > 3:
      counted = print_hit(args,hithash)
      if counted: hitcount += 1
      sys.stderr.write( "\r     queries processed:  %01.2f%%   |   hits caught:  %d" %( 100.0*count/total, hitcount ))
      query = hithash.get('query')
      querylength = hithash.get('query_length')
      numberofhits = hithash.get('numberofhits')
      hithash.clear()
      hithash['query'] = query
      hithash['query_length'] = querylength
      hithash['numberofhits'] = numberofhits
    # query
    if line.startswith('Query='):
      count += 1
      hithash['query'] = re.search('Query=\s*(\S+)',line).group(1)
      hithash['numberofhits'] = 0
    # query with no hit
    elif re.search('No hits found',line):
      if not args.has_key('evalthresh') and not args.has_key('minlength'):
        
        if args.get('delimiter') in [ 't', 'tab' ]:
          print string.join([hithash.get('query'), "no_hit_found"], '\t')
        elif args.get('delimiter') in [ ';' ',' '|' ]:
          print string.join([hithash.get('query'), "no_hit_found"], args.get('delimiter'))
        else:
          print string.join([hithash.get('query'), "no_hit_found"], ' ')
        
      hithash.clear()
      continue
    
    if args.has_key('numberofbesthits') and args.get('numberofbesthits') < hithash.get('numberofhits'):
      hithash.clear()
      continue
    
    if len(hithash) < 2: continue
    
    # query length
    if re.search("\s+\(\d+\s+letters\)",line):
      hithash['query_length'] = re.search("\s+\((\d+)\s+letters\)",line).group(1)

    # hit id and descr
    if line.startswith('>'):
      hithash['numberofhits'] += 1
      #while line != None and not re.search('Length =',line):
      #  hithash['hitdescr'] += 
      #  line = fh.readline()
      hithash['hitid'] = re.match('>(\S+)', line).group(1)
      hithash['hitdescr'] = re.match('>\S+\s+(.*)', line).group(1)
    elif hithash.has_key('hitdescr') and not hithash.has_key('sbjct_length'):
      if re.search('Length =',line):
        hithash['hitdescr'] += re.search( '(.*)\s+Length =', line ).group(1).strip()
      else:
        hithash['hitdescr'] += line.lstrip().replace('\n','')
    # subject length
    if re.search('Length =',line):
      if not hithash.has_key('sbjct_length'):
        hithash['sbjct_length'] = re.search('Length =\s{0,9}(\d+)',line).group(1)
    # hit length
    if re.search('Identities =',line):
      if not hithash.has_key('hitlength'):
        hithash['hitlength'] = re.search('Identities =\s{0,9}\d+/(\d+)',line).group(1)
      if not hithash.has_key('identities'):
        hithash['identities'] = re.search('Identities =\s{0,9}\d+/\d+\s+\((\d+)%\)',line).group(1)
    # positives
    if re.search('Positives =',line):
      if not hithash.has_key('positives'):
      hithash['positives'] = re.search('Positives =\s{0,9}\d+/\d+\s+\((\d+)%\)',line).group(1)
    # gaps
    if re.search('Gaps =',line):
      if not #TODO has key bla
      hithash['gaps'] = re.search('Gaps =\s{0,9}\d+/\d+\s+\((\d+)%\)',line).group(1)
    # score
    if re.search('Score =',line):
      hithash['score'] = re.search('Score =\s{0,9}(\S+)',line).group(1)
    # evalue
    if re.search('Expect[(\d)]* =',line):
      hithash['evalue'] = re.search('Expect[(\d)]* =\s{0,9}([0-9e.-]+)',line).group(1)
      if hithash['evalue'].count('e') > 0 and not re.match( '\d', hithash['evalue'] ):
        hithash['evalue'] = '1' + hithash['evalue']
    # frame
    if re.search('Frame =',line):
      hithash['frame'] = re.search('Frame =\s{0,9}(\S+)',line).group(1)
    # strand (BLASTN)
    if re.search('Strand =',line):
      hithash['strand'] = re.search('Strand =\s{0,9}(.*)\n',line).group(1)
    # get hit positions
    if re.search('Query:\s*\d+',line):
      if not hithash.has_key('query_startpos'): 
        hithash['query_startpos'] = re.search('Query:\s*(\d+)',line).group(1)
      hithash['query_endpos'] = re.search('(\d+)\n',line).group(1)
    if re.search('Sbjct:\s*\d+',line):
      if not hithash.has_key('sbjct_startpos'): 
        hithash['sbjct_startpos'] = re.search('Sbjct:\s*(\d+)',line).group(1)
      hithash['sbjct_endpos'] = re.search('(\d+)\n',line).group(1)      
            
  if len(hithash) > 2: 
    counted = print_hit(args,hithash)
    if counted: hitcount += 1
  sys.stderr.write( "\r     queries processed:  %01.2f%%   |   hits caught:  %d\n" %( 100.0*count/total, hitcount ))
  fh.close()
    
  
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
parse_blast_out( args )
