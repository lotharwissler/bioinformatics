#!/usr/bin/python

import os, sys 		# low level handling, such as command line stuff
import string			# string methods available
import re					# regular expressions
import getopt			# comand line argument handling
from low import *	# custom functions, written by myself
from Bio import SeqIO # biopython stuff, to parse fasta files for instance
from rpy import r
from pylab import *

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <path> -n <path>" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        path to the predicted protein sequence fasta file" )
	stdout( " -n        path to the nucleotide sequence fasta file" )
	stdout( " " )
	
	sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hf:n:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
		
	orffile, ntfile = '', ''
	for key, value in keys:
		if key == '-f': orffile = value
		if key == '-n': ntfile = value
			
	if orffile == '':
		stderr( "orf sequence data file missing." )
		show_help()
	elif not file_exists( orffile ):
		stderr( "invalid path in orffile " + orffile )
		show_help()
		
	if ntfile == '':
		stderr( "nucleotide sequence data file missing." )
		show_help()
	elif not file_exists( ntfile ):
		stderr( "invalid path in ntfile " + ntfile )
		show_help()
		
	ntfile = get_global_path( ntfile )
	orffile = get_global_path( orffile )
	return orffile, ntfile

# =============================================================================
def orf_stats( orffile, ntfile ):
	"""
	
	"""
	
	# read in all nt sequences and store them in a hash
	nthash = {}
	handle = open( ntfile )
	for seq_record in SeqIO.parse(handle, "fasta"):
		nthash[seq_record.id] = seq_record.seq.tostring()
	handle.close()
	#print "read in %s nucleotide sequences." % len(nthash)
	
	stopcodons = [ 'TAG', 'TAA', 'TGA' ]
	
	# do stats on each predicted orf (= entry in orffile)
	handle = open( orffile )
	aaseqlength = []
	aaseqstop = []
	fw = open( get_basename( orffile ) + '.withstop', 'w' )
	for seq_record in SeqIO.parse(handle, "fasta") :
		orfinfo = seq_record.description.split() # id frame ntfrom ntto
		id = seq_record.id
		aaseq = seq_record.seq.tostring()
		ntseq = nthash[id]
		
		stop = 0
		pos = int(orfinfo[-2])-1
		while (pos < int(orfinfo[-1])+3):
			codon = ntseq[ pos : pos+3 ]
			if codon in stopcodons:
				stop = 1
				break
			pos += 3
		
		print "%s\t%s\t%s" %( id, len(aaseq), stop )
		aaseqlength.append( len(aaseq) )
		aaseqstop.append( stop )
		
		if stop: fw.write( ">" + id + "\n" + aaseq + "*\n" )
			
	fw.flush()
	fw.close()
	handle.close()
	rc( 'figure', figsize=(12,5) )
	# all SNPs

	figure()
	subplot(121)
	hist( aaseqlength, fc='grey' )
	title( get_basename(orffile) + ' sequence length (aa)' )
	
	subplot(122)
	hist( aaseqstop, bins=[0,1], fc='grey' )
	title( get_basename(orffile) + ' sequences containing a stop codon' )
	
	savefig( get_basename(orffile) + '.pdf')
	

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

def main():
	"""
	"""
	orffile, ntfile = handle_arguments() 
	orf_stats( orffile, ntfile )

# =============================================================================
main()