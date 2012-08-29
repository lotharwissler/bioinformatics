#!/usr/bin/python

import os, sys 		# low level handling, such as command line stuff
from low import write_to_file
import string			# string methods available
import re					# regular expressions
import getopt			# comand line argument handling
import math				# match functions
from low import *	# custom functions, written by myself
import tempfile   # generate tmp files
from Bio import SeqIO # biopython stuff, to parse fasta files for instance

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> [-m <path> -c <path>] [-s <path>]" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        path to the fasta file containing all contig and singleton sequences" )
	stdout( " " )
	
	sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hf:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
		
	file = ''
	for key, value in keys:
		if key == '-f': file = value
			
	if file == '':
		stderr( "sequence data file missing." )
		show_help()
	elif not file_exists( value ):
		stderr( "invalid path in " + key )
		show_help()
		
	file = get_global_path( file )
	return file

# =============================================================================
def get_sequences( fastafile, number ):
	"""
	gets the first <number> of sequences within the fasta file.
	writes it to a file, returns the filename of this file.
	"""
	fh, tmpfilename = tempfile.mkstemp(dir='.')
	fw = open( tmpfilename, 'w' )
	handle = open(fastafile)
	count = 0
	for seq_record in SeqIO.parse(handle, "fasta"):
		count += 1
		if count > number: break
		fw.write( '>' + seq_record.id + '\n' + seq_record.seq.tostring() + '\n' )
	handle.close()
	fw.flush()
	fw.close()
	return tmpfilename

def generate_datasets( file ):
	
	datahash = {}
	fo = open( file )
	for line in fo:
		alphabet, counter, path = line.split()
		datahash[ alphabet+'-'+counter ] = path
	fo.close()
	
	seqhash = {}
	# nucleotide
	seqhash[ 'nt-1-1' ] = get_sequences(datahash[ 'nt-1' ], 500)
	seqhash[ 'nt-2-1' ] = get_sequences(datahash[ 'nt-2' ], 50)
	#seqhash[ 'nt-1-2' ] = get_sequences(datahash[ 'nt-1' ], 100)
	#seqhash[ 'nt-2-2' ] = get_sequences(datahash[ 'nt-2' ], 250)
	#seqhash[ 'nt-1-3' ] = get_sequences(datahash[ 'nt-1' ], 165)
	#seqhash[ 'nt-2-3' ] = get_sequences(datahash[ 'nt-2' ], 165)
	# amino acid
	#seqhash[ 'aa-1-1' ] = get_sequences(datahash[ 'aa-1' ], 100)
	#seqhash[ 'aa-2-1' ] = get_sequences(datahash[ 'aa-2' ], 100)
	#seqhash[ 'aa-1-2' ] = get_sequences(datahash[ 'aa-1' ], 150)
	#seqhash[ 'aa-2-2' ] = get_sequences(datahash[ 'aa-2' ], 150)
	#seqhash[ 'aa-1-3' ] = get_sequences(datahash[ 'aa-1' ], 200)
	#seqhash[ 'aa-2-3' ] = get_sequences(datahash[ 'aa-2' ], 200)
	#seqhash[ 'aa-1-4' ] = get_sequences(datahash[ 'aa-1' ], 300)
	#seqhash[ 'aa-2-4' ] = get_sequences(datahash[ 'aa-2' ], 300)
	
	for key, path in seqhash.iteritems():
		if key.startswith('nt'): t = 'n'
		else: t = 'p'
		os.system( "xdformat -" + t + " " + path + " &> xdformat.log")
		os.system( "formatdb -i " + path )
	
	return seqhash 

# =============================================================================
def determine_blast_program( type1, type2 ):
	"""
	"""
	if type1 == 'nt' and type2 == 'nt':
		return 'tblastx'
	elif type1 == 'aa' and type2 == 'aa':
		return 'blastp'
	elif type1 == 'aa' and type2 == 'nt':
		return 'tblastn'
	elif type1 == 'nt' and type2 == 'aa':
		return 'blastx'
	else:
		return None

# =============================================================================
def benchmark_blastall( type1, path1, type2, path2 ):
	"""
	determines the runtime of dataset1 blasted against dataset2.
	determines the type of blast to use depending on the file types (aa or nt).
	"""
	p = determine_blast_program( type1, type2 )
	starttime = time.time()
	os.system( "blastall -p " + p + " -d " + path2 + " -i " + path1 + " -o blastall.out" )
	runtime = time.time() - starttime
	print "benchmark blastall", type1, "vs", type2, "--", runtime
# =============================================================================
def benchmark_wublast( type1, path1, type2, path2 ):
	"""
	determines the runtime of dataset1 blasted against dataset2.
	determines the type of blast to use depending on the file types (aa or nt).
	wublast syntax: <program> <database> <query> [options...]
	"""
	p = determine_blast_program( type1, type2 )
	starttime = time.time()
	os.system( p + " " + path2 + " " + path1 + " &> wublast.out")
	runtime = time.time() - starttime
	sin,sout = os.popen2("grep \> -c " + path1)
	sin.close()
	s1 = sout.read().replace('\n','')
	
	sout.close()
	sin,sout = os.popen2("grep \> -c " + path2)
	sin.close()
	s2 = sout.read().replace('\n','')
	sout.close()
	
	print "benchmark wublast", s1, type1, "vs", s2, type2, "--", runtime
	
	
def xdformat( file, type ):
	os.system( "xdformat -" + type + " " + file + " &> xdformat.log")

# =============================================================================
def bench_nt_vs_aa( seqhash ):
	print "benchmark nt vs aa"
	
	ricent = 'data/rice.nt'
	riceaa = 'data/rice.aa'
	arathnt = 'data/arath.nt'
	arathaa = 'data/arath.aa'
	
	rice_nt_100 = get_sequences(ricent, 100)
	xdformat( rice_nt_100, 'n' )
	arath_nt_100 = get_sequences(arathnt, 100)
	xdformat( arath_nt_100, 'n' )
	rice_nt_300 = get_sequences(ricent, 300)
	xdformat( rice_nt_300, 'n' )
	arath_nt_300 = get_sequences(arathnt, 300)
	xdformat( arath_nt_300, 'n' )
	rice_nt_500 = get_sequences(ricent, 500)
	xdformat( rice_nt_500, 'n' )
	arath_nt_500 = get_sequences(arathnt, 500)
	xdformat( arath_nt_500, 'n' )
	
	rice_aa_100 = get_sequences(riceaa, 100)
	xdformat( rice_aa_100, 'p' )
	arath_aa_100 = get_sequences(arathaa, 100)
	xdformat( arath_aa_100, 'p' )
	rice_aa_300 = get_sequences(riceaa, 300)
	xdformat( rice_aa_300, 'p' )
	arath_aa_300 = get_sequences(arathaa, 300)
	xdformat( arath_aa_300, 'p' )
	rice_aa_500 = get_sequences(riceaa, 500)
	xdformat( rice_aa_500, 'p' )
	arath_aa_500 = get_sequences(arathaa, 500)
	xdformat( arath_aa_500, 'p' )
	
	print "---"
	print "TBLASTX"
	benchmark_wublast( 'nt', rice_nt_100, 'nt', arath_nt_100 )
	benchmark_wublast( 'nt', rice_nt_300, 'nt', arath_nt_300 )
	benchmark_wublast( 'nt', rice_nt_500, 'nt', arath_nt_500 )
	print "---"
	print "BLASTX"
	benchmark_wublast( 'nt', rice_nt_100, 'aa', arath_aa_100 )
	benchmark_wublast( 'nt', rice_nt_300, 'aa', arath_aa_300 )
	benchmark_wublast( 'nt', rice_nt_500, 'aa', arath_aa_500 )
	print "---"
	print "TBLASTN"
	benchmark_wublast( 'aa', rice_aa_100, 'nt', arath_nt_100 )
	benchmark_wublast( 'aa', rice_aa_300, 'nt', arath_nt_300 )
	benchmark_wublast( 'aa', rice_aa_500, 'nt', arath_nt_500 )
	print "---"
	print "BLASTP"
	benchmark_wublast( 'aa', rice_aa_100, 'aa', arath_aa_100 )
	benchmark_wublast( 'aa', rice_aa_300, 'aa', arath_aa_300 )
	benchmark_wublast( 'aa', rice_aa_500, 'aa', arath_aa_500 )
	print "---"
	
# =============================================================================
def bench_sizes( seqhash ):
	print "benchmark sizes"
	
	ricent = 'data/rice.nt'
	riceaa = 'data/rice.aa'
	arathnt = 'data/arath.nt'
	arathaa = 'data/arath.aa'
	
	arath_aa_200 = get_sequences(arathaa, 200)
	xdformat( arath_aa_200, 'p' )
	
	rice_aa_10 = get_sequences(riceaa, 10)
	xdformat( rice_aa_10, 'p' )
	
	rice_aa_50 = get_sequences(riceaa, 50)
	xdformat( rice_aa_50, 'p' )
	
	rice_aa_200 = get_sequences(riceaa, 200)
	xdformat( rice_aa_200, 'p' )
	
	rice_aa_300 = get_sequences(riceaa, 300)
	xdformat( rice_aa_300, 'p' )
	
	rice_aa_500 = get_sequences(riceaa, 500)
	xdformat( rice_aa_500, 'p' )
	
	print "---"
	
	benchmark_wublast( 'aa', rice_aa_10, 'aa', arath_aa_200 )
	benchmark_wublast( 'aa', rice_aa_50, 'aa', arath_aa_200 )
	benchmark_wublast( 'aa', rice_aa_200, 'aa', arath_aa_200 )
	benchmark_wublast( 'aa', rice_aa_300, 'aa', arath_aa_200 )
	benchmark_wublast( 'aa', rice_aa_500, 'aa', arath_aa_200 )
	
	print "---"
	
	benchmark_wublast( 'aa', arath_aa_200, 'aa', rice_aa_10 )
	benchmark_wublast( 'aa', arath_aa_200, 'aa', rice_aa_50 )
	benchmark_wublast( 'aa', arath_aa_200, 'aa', rice_aa_200 )
	benchmark_wublast( 'aa', arath_aa_200, 'aa', rice_aa_300 )
	benchmark_wublast( 'aa', arath_aa_200, 'aa', rice_aa_500 )
	
	print "---"
	

# =============================================================================
def bench_single_vs_multiple_files( seqhash ):
	
	def single_files( file ):
		count = 0
		filenames = {}
		handle = open(file)
		for seq_record in SeqIO.parse(handle, "fasta") :
			filenames[ file+str(count) ] = 1 
			write_to_file( file+str(count), seq_record.id + '\n' + seq_record.seq.tostring() + '\n' )
			count += 1
		handle.close()
		return filenames
			
	print "benchmark query files"
	
	ricent = 'data/rice.nt'
	riceaa = 'data/rice.aa'
	arathnt = 'data/arath.nt'
	arathaa = 'data/arath.aa'
	
	rice_aa_10 = get_sequences(riceaa, 50)
	rice_aa_50 = get_sequences(riceaa, 200)
	rice_aa_100 = get_sequences(riceaa, 500)
	arath_aa_1000 = get_sequences(arathaa, 1000)
	xdformat( arath_aa_1000, 'p' )
	
	print "---"
	
	# split the files
	p = 'blastp'

	
	starttime = time.time()
	filenames = single_files( rice_aa_10 )
	for file in filenames.keys():
		os.system( p + " " + arath_aa_1000 + " " + file + " &> wublast.out")
		sys.stdout.write('.')
	runtime = time.time() - starttime
	sys.stdout.write('\n')
	print "benchmark wublast", str(len(filenames.keys())), "--", runtime
	for file in filenames.keys(): os.unlink(file)
	
	starttime = time.time()
	filenames = single_files( rice_aa_50 )
	for file in filenames.keys():
		os.system( p + " " + arath_aa_1000 + " " + file + " &> wublast.out")
		sys.stdout.write('.')
	runtime = time.time() - starttime
	sys.stdout.write('\n')
	print "benchmark wublast", str(len(filenames.keys())), "--", runtime
	for file in filenames.keys(): os.unlink(file)
	
	starttime = time.time()
	filenames = single_files( rice_aa_100 )
	for file in filenames.keys():
		os.system( p + " " + arath_aa_1000 + " " + file + " &> wublast.out")
		sys.stdout.write('.')
	runtime = time.time() - starttime
	sys.stdout.write('\n')
	print "benchmark wublast", str(len(filenames.keys())), "--", runtime
	for file in filenames.keys(): os.unlink(file)
	
	print "---"
	benchmark_wublast( 'aa', rice_aa_10, 'aa', arath_aa_1000 )
	benchmark_wublast( 'aa', rice_aa_50, 'aa', arath_aa_1000 )
	benchmark_wublast( 'aa', rice_aa_100, 'aa', arath_aa_1000 )

# =============================================================================
def remove_tmpfiles( seqhash ):
	for key, value in seqhash.iteritems():
		os.system( "rm " + value + "*" )

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

def main():
	"""
	"""
	file = handle_arguments() 
	seqhash = generate_datasets( file )
	#bench_nt_vs_aa( seqhash )
	bench_sizes( seqhash )
	#bench_single_vs_multiple_files( seqhash )
	remove_tmpfiles( seqhash ) 

# =============================================================================
main()