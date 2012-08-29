#!/usr/bin/python

import os, sys 		# low level handling, such as command line stuff
import string			# string methods available
import re					# regular expressions
import getopt			# comand line argument handling
import anydbm					# index databases (file hash)
from low import *	# custom functions, written by myself
from Bio import SeqIO # biopython stuff, to parse fasta files for instance

class GOterm:
	def __init__(self):
		self.goid = ""
		self.name = ""
		self.namespace = ""
		self.definition = ""
		self.synonym = ""
		self.is_a = []
		self.alt_id = []
		
	def tostring(self):
		S = ""
		S += "id: " + self.goid
		S += "\tname: " + self.name
		S += "\tnamespace: " + self.namespace
		S += "\tdef: " + self.definition
		S += "\tsynonym: " + self.synonym
		S += "\talt_id: " + string.join(self.alt_id,',')
		S += "\tis_a: " + string.join(self.is_a,',')
		return S

class OBOParser:
	
	def __init__(self,file):
		self.file = file
		self.filehandle = open( file )
		self.alt_ids = {}
		
	def next(self):
		ngo = None
		for line in self.filehandle:
			if line.startswith('[Term]'):	ngo = GOterm()
			if ngo == None: continue
			if line.startswith('id:'): ngo.goid = re.search( '^id:\s+(.*)$', line ).group(1)
			if line.startswith('name:'): ngo.name = re.search( '^name:\s+(.*)$', line ).group(1)
			if line.startswith('namespace:'): ngo.namespace = re.search( '^namespace:\s+(.*)$', line ).group(1)
			if line.startswith('def:'): ngo.definition = re.search( '^def:\s+"(.*)"', line ).group(1)
			if line.startswith('synonym:'): ngo.synonym = re.search( '^synonym:\s+"(.*)"', line ).group(1)
			if line.startswith('alt_id:'): self.alt_ids[ re.search( '^alt_id:\s+(.*)$', line ).group(1) ] = ngo.goid
			if line.startswith('is_a:'): ngo.is_a.append( re.search( '^is_a:\s+(.*)$', line ).group(1) )
			if re.match( '$', line): break
		return ngo
		
	def get_alt_ids(self):
		return self.alt_ids
	
	def close(self):
		self.filehandle.close()
	
	
		

# =============================================================================	
def show_help( ):
	""" displays the program parameter list and usage information """
	stdout( "usage: " + sys.argv[0] + " -f <oath> -o <path> [ -e ]" )
	stdout( " " )
	stdout( " option    description" )
	stdout( " -h        help (this text here)" )
	stdout( " -f        parsed gene ontology blast file" )
	stdout( " -o        OBO file" )
	stdout( " -e        evalue threshold. default: 10.0" )
	stdout( " " )
	sys.exit(1)

# =============================================================================
def handle_arguments():
	""" verifies the presence of all necessary arguments and returns the data dir """
	if len ( sys.argv ) == 1:
		stderr( "no arguments provided." )
		show_help()	
	
	try: # check for the right arguments
		keys, values = getopt.getopt( sys.argv[1:], "hf:o:e:" )
	except getopt.GetoptError:
		stderr( "invalid arguments provided." )
		show_help()
	
	args = {}
	for key, value in keys:
		if key == '-f':	args['file'] = value
		if key == '-o':	args['OBO'] = value
		if key == '-e':	args['evalue'] = float(value)
				
	if not args.has_key('file'):
		stderr( "GO file missing." )
		show_help()
	if not file_exists( args.get('file') ):
		stderr( "GO file does not exist." )
		show_help()

	if not args.has_key('OBO'):
		stderr( "OBO file missing." )
		show_help()
	if not file_exists( args.get('OBO') ):
		stderr( "OBO file does not exist." )
		show_help()
		
	if not args.has_key('evalue'):
		args['evalue'] = 10.0
		
	return args


# =============================================================================
def get_obo_hash( file ):
	obohash = {}
	OP = OBOParser( file )
	gotermcount = 0
	while (1):
		goterm = OP.next()
		if goterm == None: break
		obohash[ goterm.goid ] = goterm
		gotermcount += 1
		sys.stderr.write( "\r     processing OBO file ...   |   goterms caught:  %d" %(gotermcount) )
			
	alt_ids = OP.get_alt_ids()
	for altid, goid in alt_ids.iteritems():
		obohash[ altid ] = obohash.get(goid)
	
	return obohash
	OP.close()
	
# =============================================================================
# =============================================================================
def main( args ):
	
	obohash = get_obo_hash( args.get('OBO') )
	
	oldquery = ''
	querygoterms = {}
	regex = re.compile('GO:\d+')
	entrycount = 0
	gotermcount = 0
	sout, serr = catch_bash_cmd_output("cat %s | wc -l" %args.get('file'))
	totalentries = int(sout)
	
	fo = open( args.get('file') )
	for line in fo:
		line = line.replace('\n','')
		columns = line.split('\t')
		entrycount += 1
		sys.stderr.write( "\r     entries processed:  %01.2f%%   |   goterms caught:  %d" %( 100.0*entrycount/totalentries, gotermcount ))
		# not hit found: skip
		if columns[1] == 'no_hit_found': continue
		# else
		qid, hitid, evalue, descr = columns
		if oldquery != qid:
			oldquery = qid
			querygoterms.clear()
			print ">%s" %qid
		assoc_goterms = re.findall( regex, descr )
		for e in assoc_goterms:
			if querygoterms.has_key(e): continue
			querygoterms[e] = 1
			if float(evalue) < args.get('evalue'):
				List = [e]
				if obohash.has_key(e):
					List.append( obohash.get(e).namespace )
					List.append( obohash.get(e).name )
					#List.append( obohash.get(e).namespace )
				else:
					stderr( "GO:id not found in the OBO hash: \"%s\"" %e )
				List.append(evalue)
				List.append(hitid)
				print string.join( List, '\t' )
				gotermcount += 1
	fo.close()
	sys.stderr.write( "\r     entries processed:  %01.2f%%   |   goterms caught:  %d\n" %( 100.0*entrycount/totalentries, gotermcount ))
	
	
# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )