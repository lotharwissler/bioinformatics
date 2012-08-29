#!/usr/bin/python

import os, sys 		# low level handling, such as command line stuff
import string			# string methods available
import re					# regular expressions
import getopt			# comand line argument handling
import math				# match functions
from low import *	# custom functions, written by myself



# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -b <path> [-f <path>]" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        path to the fasta file containing the record ids." )
	stdout( " -b        path to the blast.best-hit file of swiss-prot" )
	stdout( " " )
	
	sys.exit(1)


# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hb:f:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
		
	blastbesthitfile, recordfile = '', ''
	for key, value in keys:
		if key == '-b':
			if not file_exists( value ):
				stderr( "invalid path in " + key )
				show_help()
			else:
				blastbesthitfile = value
		
		if key == '-f':
			if not file_exists( value ):
				stderr( "invalid path in " + key )
				show_help()
			else:
				recordfile = value
		
	if blastbesthitfile == '':
		stderr( "blast.best-hit file missing." )
		show_help()
	elif not file_exists( blastbesthitfile ):
		stderr( "blast.best-hit file does not exist." )
		show_help()
		
	if recordfile == '':
		stderr( "recordfile missing." )
		show_help()
	elif not file_exists( recordfile ):
		stderr( "recordfile does not exist." )
		show_help()
		
	return blastbesthitfile, recordfile


# =============================================================================
def parse_best_blast_hits( blastbesthitfile, recordfile ):
	""" """
	
	records = []
	fo = open( recordfile, 'r' )
	for line in fo:
		records.append(line.strip().replace('\n',''))
	
	fo = open( blastbesthitfile, 'r' )
	for line in fo:
		columns = line.split()
		if columns[0] in records:
			print columns[0]
			print "   hit   :", string.join(columns[10:], ' ')[1:]
			print "   evalue:", columns[4], "\n"
		
	fo.close()
	
	

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

blastbesthitfile, recordfile = handle_arguments()
parse_best_blast_hits( blastbesthitfile, recordfile )
