#!/usr/bin/python

import os, sys 		# low level handling, such as command line stuff
import string			# string methods available
import re					# regular expressions
from low import *	# custom functions, written by myself


# =============================================================================
def get_translatedfasta_from_gb( file ):
	"""
	"""
	
	def write_output( source, hash ):
		L = [ ">",hash.get('protein_id'),"|",hash.get('db_xref')]
		if hash.has_key('product'): L.append("|"+hash.get('product'))
		L.append(" ("+source+")")
		print string.join(L,'')
		print hash.get('translation')
	
	fo = open(file)
	# read general infos
	source = ''
	for line in fo:
		if re.search('FEATURES',line): break
		if re.match('SOURCE',line):
			source = re.search('SOURCE\s+(.*)\n',line).group(1)
		
	# read gene infos
	hash = {}
	hit = 0
	for line in fo:
		if not re.match('                     ',line):
			if len(hash) > 0:	write_output( source, hash )
			hash = {} 
			hit = 0
		if re.match('     CDS',line): hit = 1
		
		if hit:
			# catch everything except translation sequence
			if re.search('/(\S+)=".*"',line):
				hash[re.search('/(\S+)=".*"',line).group(1)] = re.search('/\S+="(.*)"',line).group(1)
			# catch translation sequence
			if re.search('/translation=',line):
				hash['translation'] = re.search('/translation="(.*)\n',line).group(1)
			elif hash.has_key('translation'): hash['translation'] += re.search("([a-zA-Z]+)",line).group(1)
	if len(hash) > 0: write_output( source, hash )
	fo.close()

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

if len( sys.argv ) == 1:
	print "no arguments provided. you need to specify the gb file(s) to parse."
	sys.exit(1)

for file in sys.argv[1:]:
	if not file_exists(file):
		print "gb file not found (or is a dir):", file
		continue
	get_translatedfasta_from_gb( file ) 