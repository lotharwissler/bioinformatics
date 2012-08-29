#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself
import anydbm					# index databases (file hash)
from Bio import SeqIO # biopython stuff, to parse fasta files for instance

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> -s <path> -d <path> -o <path> [-i]" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        in cluster file" )
	stdout( " -s        swiss-prot database folder" )
	stdout( " -d        genome datasets folder" )	
	stdout( " -i        reindex databases" )
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
		keys, values = getopt.getopt( sys.argv[1:], "hf:o:d:s:i" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {}
	for key, value in keys:
		if key == '-f': args['in'] = value
		if key == '-o':	args['out'] = value
		if key == '-d':	args['datasets'] = value
		if key == '-s':	args['swissprot'] = value
		if key == '-i':	args['indexdb'] = 1
				
	if not args.has_key('in'):
		stderr( "in file missing." )
		show_help()
	if not file_exists( args.get('in') ):
		stderr( "in file does not exist." )
		show_help()
		
	if not args.has_key('datasets'):
		stderr( "datasets folder missing." )
		show_help()
	if not dir_exists( args.get('datasets') ):
		stderr( "datasets folder does not exist." )
		show_help()
		
	if not args.has_key('swissprot'):
		stderr( "swissprot folder missing." )
		show_help()
	if not dir_exists( args.get('swissprot') ):
		stderr( "swissprot folder does not exist." )
		show_help()
		
		
	if not args.has_key('out'):
		stderr( "out folder missing." )
		show_help()
	
	if not dir_exists( args.get('out') ):
		os.mkdir( args.get('out') )
	
	if not args['out'].endswith('/'): args['out'] += '/'
	if not args['swissprot'].endswith('/'): args['swissprot'] += '/'
	if not args['datasets'].endswith('/'): args['datasets'] += '/'
	
	return args

# =============================================================================
def get_cluster_ids( file ):
	hash = {}
	count = 0
	fo = open( file, 'r' )
	for line in fo:
		count += 1
		hash[count] = line.split()
	fo.close()
	return hash


# =============================================================================
def get_fasta_( id, args, genomes ):
	if id[:2] in genomes:
		db = args.get('datasets') + id[:2] + '.aa'
	else: db = args.get('swissprot') + 'uniprot_sprot.fasta'
	print "lookup | db:", db, "\t->", id
	
	if not file_exists( db + '.xpi' ):
		os.system( "xdformat -p -I %s &> /dev/null" % db )
		
	out = os.popen( "xdget -p %s %s" %( db, id ) )
	seq = out.read()
	out.close()
	return seq
	

# =============================================================================
def index_databases( dbs, args ):
	ext = os.path.splitext(dbs[0])[1]
	DBM_name = '/data/l_wiss01/database/all-fasta-records.dbm' + ext
	if file_exists(DBM_name) and not args.has_key('indexdb'):
		return DBM_name
	print "creating DBM:",  DBM_name
	DBM = anydbm.open( DBM_name, 'c' )
	for db in dbs:
		print "-> adding db:", db
		handle = open(db)
		for seq_record in SeqIO.parse(handle, "fasta"):
			DBM[ seq_record.id ] = seq_record.seq.tostring()
		handle.close()
	DBM.close()
	print "DONE. indexed database:", DBM_name
	return DBM_name
	
# =============================================================================
# =============================================================================
def main( args ):
	genomes = ['PO', 'ZO', 'At', 'Os', 'Mt', 'Pt', 'Lj' ]
	# database names
	aa_dbs = []
	for g in genomes:
		aa_dbs.append( args.get('datasets') + g + '.aa' )
	nt_dbs = []
	for g in genomes:
		nt_dbs.append( args.get('datasets') + g + '.nt' )
	#dbs.append( args.get('swissprot') + + 'uniprot_sprot.fasta' )
	
	# index databases
	aa_dbmname = index_databases( aa_dbs, args )
	nt_dbmname = index_databases( nt_dbs, args )
	
	clusterhash = get_cluster_ids( args.get('in') )
	
	aa_db = anydbm.open(aa_dbmname, "r")
	nt_db = anydbm.open(nt_dbmname, "r")
	
	for i, idlist in clusterhash.iteritems():
		fwaa = open( args.get('out') + 'cluster' + add_leading_zeroes(i,3) + '.aa', 'w' )
		fwnt = open( args.get('out') + 'cluster' + add_leading_zeroes(i,3) + '.nt', 'w' )
		for id in idlist[1:]:
			if not aa_db.has_key(id) or not nt_db.has_key(id):
				stderr( "cluster %s | id %s not in both datasets | skipped." %(i, id) )
				continue
			fwaa.write( ">" + id + "\n" + aa_db[ id ] + "\n" )
			fwnt.write( ">" + id + "\n" + nt_db[ id ] + "\n" )
		fwaa.flush()
		fwaa.close()
		fwnt.flush()
		fwnt.close()
		
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )