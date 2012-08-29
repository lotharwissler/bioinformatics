#!/usr/bin/python
import os, sys
import rpy
import string
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself



# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -a <path> -t <path> -m <N> -n <namespaces>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -a        annotation file in topGO table format" )
  stdout( " -t        test ids file" )
  stdout( " -e        use ELIM instead of WEIGHT algorithm" )
  stdout( " -n        list of namespaces to test. default: \"BP,CC,MF\"" )
  stdout( " -m        minimum number of genes per GO term for the GO term to be tested. default: 1" )
  stdout( " -o        test for over-representation" )
  stdout( " -u        test for under-representation" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "ha:t:n:m:oue" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'namespaces':["BP", "CC", "MF"], 'min':1, 'over':False, 'under':False, 'algorithm':'weight01'}
  for key, value in keys:
    if key == '-a': args['annot'] = value
    if key == '-t': args['testset'] = value
    if key == '-m': args['min'] = int(value)
    if key == '-n': args['namespaces'] = value.split(",")
    if key == '-o': args['over'] = True
    if key == '-u': args['under'] = True
    if key == '-e': args['algorithm'] = 'elim'
    
  if not args.has_key('annot'):
    stderr( "annot file argument missing." )
    show_help()
  elif not file_exists( args.get('annot') ):
    stderr( "annot file does not exist." )
    show_help()
  
  if not args.has_key('testset'):
    stderr( "testset file argument missing." )
    show_help()
  elif not file_exists( args.get('testset') ):
    stderr( "testset file does not exist." )
    show_help()

  return args


# =============================================================================
def init_R(under=False):
  R = rpy.r
  R('sink("/dev/null")')
  try:
    R('invisible(capture.output(library("topGO")))')
  except:
    try: 
      R.source("http://bioconductor.org/biocLite.R")
      R.biocLite('topGO')
      R.library('topGO')
    except:
      print "Problem importing R libraries."
      sys.exit()
  
  if under:
    R('if(!isGeneric("GOFisherTestUnder")) setGeneric("GOFisherTestUnder", function(object) standardGeneric("GOFisherTestUnder"))')
    R('setMethod("GOFisherTestUnder", "classicCount", function(object) { contMat <- contTable(object); if(all(contMat == 0)) p.value <- 1 else p.value <- fisher.test(contMat, alternative = "less")$p.value; return(p.value) })')
  return R

# =============================================================================
def test_for_overrepresentation(R, args):
  significant = []
  tmp = R('GOmap = readMappings(file = "' + args['annot'] + '")')
  tmp = R('refset  = names(GOmap)')
  tmp = R('testset = scan(file="' + args['testset'] + '", what=character())')
  tmp = R('genes_of_interest = factor(as.integer(refset %in% testset))')
  tmp = R('names(genes_of_interest) <- refset')
  for ontology in args['namespaces']:
    tmp = R('tgData = new("topGOdata", ontology = "' + ontology + '", allGenes = genes_of_interest, nodeSize = ' + str(args['min']) + ', annot = annFUN.gene2GO, gene2GO = GOmap)')
    pvalueHash = R('score(runTest(tgData, algorithm="%s", statistic="fisher"))' %(args['algorithm']))
    keys, pvalues = [], []
    for key, p in pvalueHash.iteritems():
      keys.append(key)
      pvalues.append(p)
    tmp = R.assign('pvalues',pvalues)
    padjusted = R('p.adjust(pvalues, method="fdr")')
    for i in range(len(keys)):
      goterm = keys[i]
      p = pvalues[i]
      fdr = padjusted[i]
      if p > 0.05: continue
      significant.append(["O", ontology, goterm, str(p), str(fdr)])
  return significant

# =============================================================================
def test_for_underrepresentation(R, args):
  significant = []
  tmp = R('GOmap = readMappings(file = "' + args['annot'] + '")')
  tmp = R('refset  = names(GOmap)')
  tmp = R('testset = scan(file="' + args['testset'] + '", what=character())')
  tmp = R('genes_of_interest = factor(as.integer(refset %in% testset))')
  tmp = R('names(genes_of_interest) <- refset')
  for ontology in args['namespaces']:
    tmp = R('tgData = new("topGOdata", ontology = "' + ontology + '", allGenes = genes_of_interest, nodeSize = ' + str(args['min']) + ', annot = annFUN.gene2GO, gene2GO = GOmap)')
    tmp = R('test.stat <- new("weightCount", testStatistic = GOFisherTestUnder, name ="Fisher test underrepresentation")')
    pvalueHash = R('score(getSigGroups(tgData, test.stat))')
    keys, pvalues = [], []
    for key, p in pvalueHash.iteritems():
      keys.append(key)
      pvalues.append(p)
    tmp = R.assign('pvalues',pvalues)
    padjusted = R('p.adjust(pvalues, method="fdr")')
    for i in range(len(keys)):
      goterm = keys[i]
      p = pvalues[i]
      fdr = padjusted[i]
      if p > 0.05: continue
      significant.append(["U", ontology, goterm, str(p), str(fdr)])
  return significant


# =============================================================================
def main(args):
  
  R = init_R(args['under'])
  fw = open(args['testset'] + ".ORA", "w")
  if args['over']: 
    results = test_for_overrepresentation(R, args)
    for r in results: fw.write(string.join(r, "\t") + "\n")
  if args['under']: 
    results = test_for_underrepresentation(R, args)
    for r in results: fw.write(string.join(r, "\t") + "\n")
  fw.close()



args = handle_arguments()
main( args )
