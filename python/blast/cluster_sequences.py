#!/usr/bin/python

import os, sys 		# low level handling, such as command line stuff
import string			# string methods available
import re					# regular expressions
import getopt			# comand line argument handling
from low import *	# custom functions, written by myself

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> -e <1e-n> -i <75> -p <90> -o <path>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        in sequence file" )
	stdout( " -e        maximum evalue" )
	stdout( " -i        minimum percent identity" )
	stdout( " -p        minimum percent positives" )
	stdout( " -o        out folder" )
	stdout( " " )
	
	sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hf:o:e:i:p:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {}
	for key, value in keys:
		if key == '-f': args['in'] = value
		if key == '-o':	args['out'] = value
		if key == '-i':	args['identities'] = value
		if key == '-p':	args['positives'] = value
		if key == '-e':	args['evalue'] = value
				
	if not args.has_key('in'):
		stderr( "in file missing." )
		show_help()
	if not file_exists( args.get('in') ):
		stderr( "in file does not exist." )
		show_help()
		
	if not args.has_key('out'):
		stderr( "out folder missing." )
		show_help()
	
	if not dir_exists( args.get('out') ):
		os.mkdir( args.get('out') )
	
	if not args['out'].endswith('/'): args['out'] += '/'
	
	return args

# =============================================================================
def blast( query, subject, overwrite=0 ):
	#infomsg( "BLAST %s vs. %s" %(query, subject) )
	blastout = "blast-out." + get_basename(query).replace('red_','') + "_" + get_basename(subject)
	# formatdb
	#os.system( "formatdb -i %s" %(subject) )
	# blast
	if not file_exists( blastout ) or overwrite:
		os.system( "blastall -p blastp -i %s -d %s -o %s" %( query, subject, blastout ) )
	#print "blastall -p blastp -i %s -d %s -o %s" %( query, subject, blastout )
	return blastout
	
# =============================================================================
def parse_blastout( file, args ):
	parseout = "parsed" + os.path.splitext( file ) [1]
	programmcall = "parse_blast_out2.py -f " + file
	if args.has_key( 'evalue' ): programmcall += " -e " + args.get('evalue')
	if args.has_key( 'identities' ): programmcall += " -i " + args.get('identities')
	if args.has_key( 'positives' ): programmcall += " -p " + args.get('positives')
	#if not file_exists( parseout ):
	os.system( programmcall + " > " + parseout )

	# make list non-redundant: remove self hits
	nonred_out = args.get('out') + 'blast-out.parsed-nr'
	fo = open(parseout)
	fw = open( nonred_out, 'w' )
	for line in fo:
		qid, hid = line.split()[0:2]
		if qid == hid: continue
		fw.write( line )
	fw.flush()
	fw.close()
	fo.close()
	return nonred_out


# =============================================================================
def assemble_clusters( parsedfile, args ):
	hash = {}
	fo = open( parsedfile, 'r' )
	for line in fo:
		qid, hid = line.split()[0:2]
		#print qid, "\t", hid
		if not hash.has_key(qid) and not hash.has_key(hid):
			hash[ qid ] = [hid]
		else:
			if qid == hid: 
				#print "skipped"
				continue
			if hash.has_key(qid) and not hid in hash.get(qid):
				hash[ qid ].append( hid )
			elif hash.has_key(hid) and not qid in hash.get(hid):	
				hash[ hid ].append( qid )
			#else:
			#	print "skipped"
		#else:
		#	stderr( "strange hash behavior! %s %s" %( qid, hid ) )
	fo.close()
	clusterout = args.get('out') + 'clusters.ids'
	fw = open( clusterout , 'w' )
	for id, listofids in hash.iteritems():
		fw.write( id + ' ' + string.join( listofids, ' ' ) + '\n' )
	fw.flush()
	fw.close()
	return clusterout
	

# =============================================================================
def annotate_clusters( clusterout, args ):
	
	def get_sequence( id ):
		seqfile = args.get('out')+'seqfile.tmp'
		os.system( "xdget -p ~/workspace/EST/rsd/datasets/aa/%s %s > %s" %( id[0:2]+'.aa', id, seqfile ) )
		return seqfile
	
	def blast_sp( seqfile ):
		out = blast( seqfile, '~/workspace/EST/blast-databases/swissprot-all', 1 )
		return out
	
	def get_descr( blastout ):
		descr = ''
		fo = open( blastout, 'r' )
		for line in fo:
			#print line
			if descr == '':
				if not line[:1] == '>': continue
				else:	descr += line[1:].replace('\n','')
			else:
				if re.search( 'Length =', line):
					descr += line[:line.index('Length =')]
					break
				else:
					descr += line.replace('\n','')
		fo.close()
		descr = re.sub( '\s{2,99}', ' ' , descr )
		return descr
	
	fo = open( clusterout, 'r' )
	for line in fo:
		ids = line.split()
		descr = []
		for id in ids:
			seq = get_sequence( id )
			blastout = blast_sp( seq )
			descr.append( get_descr( blastout ) )
		print "cluster:", string.join(ids, ' ')
		for d in descr:
			print "descr:", d
		print " -" *50
		#sys.exit(1212)
	
	fo.close()

# =============================================================================
# =============================================================================
def main( args ):
	
	blastout = blast( args.get('in'), args.get('in') )
	parsedfile = parse_blastout( blastout, args )
	clusterout = assemble_clusters( parsedfile, args )
	annotate_clusters( clusterout, args )
	
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )