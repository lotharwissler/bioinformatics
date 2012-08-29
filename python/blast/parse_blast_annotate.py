#!/usr/bin/python

import os, sys 		# low level handling, such as command line stuff
import string			# string methods available
import re					# regular expressions
import getopt			# comand line argument handling
from low import *	# custom functions, written by myself

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> [-b <path> -o <path>]" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        blast.out file to be parsed" )
	stdout( " -n        number of best hits to parse" )
	stdout( " -e        minimum evalue of a hit to be parsed" )
	stdout( " -i        minimum identity (in %)" )
	stdout( " -l        minimum length of a hit to be parsed" )
	stdout( " -d        delimiter" )
	stdout( "           default is the blast.out base name plus .best-hit" )
	stdout( " " )
	
	sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hb:f:n:e:l:i:p:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {}
	args['delimiter'] = '\t'
	for key, value in keys:
		if key == '-f':
			if not file_exists( value ):
				stderr( "invalid path in " + key )
				show_help()
			else:
				args['file'] = value
		if key == '-e':	args['evalthresh'] = float(value)
		if key == '-l':	args['minlength'] = int(value)
		if key == '-i':	args['minident'] = int(value)
		if key == '-p':	args['minpos'] = int(value)
		if key == '-n':	args['numberofbesthits'] = int(value)
		if key == '-d':	args['delimiter'] = value
				
	if not args.has_key('file'):
		stderr( "blast out file missing." )
		show_help()
	
	return args

# =============================================================================
def print_hit(args,hithash):
	"""
	"""
	if args.has_key('evalthresh'):
		if args.get('evalthresh') < float(hithash.get('evalue')):	return
	if args.has_key('minlength'):
		if args.get('minlength') > int(hithash.get('hitlength')): return
	if args.has_key('minident'):
		if args.get('minident') > int(hithash.get('identities')): return
	if args.has_key('minpos'):
		if args.get('minpos') > int(hithash.get('positives')): return
		
	hithash['hitdescr'] = re.sub( '\s{2,99}', ' ' , hithash['hitdescr'] )
		
	L = []
	L.append(hithash.get('query'))
	L.append(hithash.get('evalue'))
	#splits = hithash.get('hitdescr').split()
	#id = splits[0]
	#descr, species = string.join( splits[1:], ' ').split(' - ')
	#L.append( id )
	#L.append( descr )
	#L.append( species )
	L.append( hithash.get('hitdescr') )
	
	print string.join(L, args.get('delimiter') )
	

# =============================================================================
def parse_blast_out( args ):
	#print "# blast.out file:", args.get('file')
	#print "# numberofbesthits:", args.get('numberofbesthits')
	#print "# max.evalue:", args.get('evalthresh')
	#print "# min.length:", args.get('minlength')
	#print "# fields: query, hitid, score, evalue, query_startpos, query_endpos, sbjct_startpos, sbjct_endpos, hitlength, length, identities, positives, frame_or_strand"
	
	hithash = {}
	currenthits = 0
	fh = open( args.get('file') )
	for line in fh:
		# new hit entry
		if (line.startswith('Query=') or line.startswith('>')) and len(hithash) > 1:
			print_hit(args,hithash)
			query = hithash.get('query')
			hithash.clear()
			hithash['query'] = query
			currenthits += 1
		# query
		if line.startswith('Query='):
			hithash['query'] = re.search('Query=\s*(\S+)',line).group(1)
			currenthits = 0
		# query with no hit
		elif re.search('No hits found',line):
			if not args.has_key('evalthresh') and not args.has_key('minlength'):
				print hithash.get('query') + args.get('delimiter') + "no_hit_found"
			hithash.clear()
			currenthits = 0
			continue
		
		if args.has_key('numberofbesthits') and not args.get('numberofbesthits') > currenthits:
			hithash.clear()
			continue
		
		if len(hithash) < 1: continue
		
		# hit id and descr
		if not hithash.has_key('hitdescr'):
			if line[:1] == '>':
				hithash['hitdescr'] = line[1:].replace('\n','')
				hithash['hitid'] = re.match('(\w)',hithash['hitdescr']).group(1)
		else:
			if not hithash.has_key('length'):
				if re.search( 'Length =', line):
					hithash['hitdescr'] += line[:line.index('Length =')]
				else:
					hithash['hitdescr'] += line.replace('\n','')
		
		# subject length
		if re.search('Length =',line):
			hithash['length'] = re.search('Length =\s{0,9}(\d+)',line).group(1)
		# hit length
		if re.search('Identities =',line):
			hithash['hitlength'] = re.search('Identities =\s{0,9}\d+/(\d+)',line).group(1)
		# identities
		if re.search('Identities =',line):
			hithash['identities'] = re.search('Identities =\s{0,9}\d+/\d+\s+\((\d+)%\)',line).group(1)
		# positives
		if re.search('Positives =',line):
			hithash['positives'] = re.search('Positives =\s{0,9}\d+/\d+\s+\((\d+)%\)',line).group(1)
		# gaps
		if re.search('Gaps =',line):
			hithash['gaps'] = re.search('Gaps =\s{0,9}\d+/\d+\s+\((\d+)%\)',line).group(1)
		# score
		if re.search('Score =',line):
			hithash['score'] = re.search('Score =\s{0,9}(\S+)',line).group(1)
		# evalue
		if re.search('Expect[(2)]* =',line):
			hithash['evalue'] = re.search('Expect[(2)]* =\s{0,9}([0-9e.-]+)',line).group(1)
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
						
	if len(hithash) > 1: print_hit(args,hithash)
	fh.close()
		
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
parse_blast_out( args )