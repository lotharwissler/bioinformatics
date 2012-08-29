#!/usr/bin/python

import os, sys 		# low level handling, such as command line stuff
import string			# string methods available
import re					# regular expressions
import getopt			# comand line argument handling
from low import *	# custom functions, written by myself

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> -e <1e-10> -i <70> -p <85> -o <path>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        in homolog source file" )
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
def get_src_files( file ):
	return read_from_file( file ).splitlines()

# =============================================================================
def blast( query, subject ):
	infomsg( "BLAST %s vs. %s" %(query, subject) )
	blastout = "blast-out." + get_basename(query).replace('red_','') + "_" + get_basename(subject)
	# formatdb
	if file_exists( subject) and not file_exists( subject + '.pin' ):
		os.system( "formatdb -i %s" %(subject) )
	# blast
	if not file_exists( blastout ):
		os.system( "blastall -p blastp -i %s -d %s -o %s" %( query, subject, blastout ) )
	return blastout
	
# =============================================================================
def parse_blastout( file, args ):
	parseout = "parsed" + os.path.splitext( file ) [1]
	programmcall = "parse_blast_out2.py -f " + file
	if args.has_key( 'evalue' ): programmcall += " -e " + args.get('evalue')
	if args.has_key( 'identities' ): programmcall += " -i " + args.get('identities')
	if args.has_key( 'positives' ): programmcall += " -p " + args.get('positives')
	if not file_exists( parseout ):
		os.system( programmcall + " > " + parseout )
	hash = {}
	fo = open(parseout)
	for line in fo:
		hash[ line.split()[1] ] = 1
	fo.close()
	return parseout, hash

# =============================================================================
def get_orthologs( hashes, id, orthologlist ):
	if not hashes.has_key( id[0:2]): return orthologlist
	searchhash = hashes.get(id[0:2])
	if not searchhash.has_key(id): return orthologlist
	orths = searchhash.get(id)
	for o in orths:
		if o in orthologlist: continue
		orthologlist.append(o)
		get_orthologs( hashes, o, orthologlist )
		
	return orthologlist

# =============================================================================
def integrate_all_homologs( files, args ):
	# gather all data into a hash per pair
	hashes = {}
	for file in files:
		hash = {}
		fo = open( file, 'r')
		for line in fo:
			qid, sid = line.split()[0:2]
			if not hash.has_key( qid ): hash[ qid ] = [sid]
			else: hash[ qid ].append( sid )
		fo.close()
		hashes[ qid[0:2] ] = hash
	
	# iterate all hashes, integrate homologous relationship into a single big hash
	# take the first hash (query ids) as the reference keyset
	fw = open( args.get('out') + 'homologous-clusters.txt', 'w' )
	Homologs = {}
	for searchid in hashes.get('PO').keys():
		orthologlist = get_orthologs( hashes, searchid, [] )
		Homologs[ searchid ] = orthologlist
		fw.write( str(len(orthologlist)+1)+' '+ searchid + ' ' + string.join(orthologlist, ' ') + '\n' )
	
	fw.flush()
	fw.close()
	return Homologs
			

# =============================================================================
# =============================================================================
def main( args ):
	
	sourcefiles = get_src_files( args.get('in') )
	hashofhos = None
	parsedfiles = []
	for i in range( len(sourcefiles)-1 ):
		
		queryfile = sourcefiles[i]
		subjectfile = sourcefiles[i+1]
		
		if hashofhos:
			idfile = 'keepids.tmp'
			fw = open( idfile, 'w' )
			for id in hashofhos.keys():	fw.write( id + '\n')			
			fw.flush()
			fw.close()
			outfile = 'red_' + get_basename(queryfile) + '.aa'
			os.system( 'reduce_fasta_file.py -f %s -i %s -o %s' %(queryfile,idfile,outfile) )
			queryfile = outfile
		
		blastout = blast( queryfile, subjectfile )
		parsedfile, hashofhos = parse_blastout( blastout, args )
		parsedfiles.append( parsedfile )
		infomsg( "hits: %s" %len(hashofhos) )
		
	Homologs = integrate_all_homologs( parsedfiles, args )
	
	# stats
	no = []
	for sid, orthlist in Homologs.iteritems():
		n = len(orthlist) + 1
#		infomsg( str(n) )
		no.append(n)
	
	from rpy import r
	outfile = 'hist_size_homol_sets.pdf'
	title = 'Size of Homologous Sets'
	x = 'number of homologs'
	y = 'frequency'
	r.pdf( outfile )
	r.hist(no, xlab=x, ylab=y, main=title, col='grey', breaks=max(no))
	r.dev_off()
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )