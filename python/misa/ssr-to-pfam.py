#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import getopt      # comand line argument handling
from collections import defaultdict
from low import *  # custom functions, written by myself
from misa import MisaSSR
from pfam import PfamDomain
import rpy2.robjects as robjects
R = robjects.r

# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  print >> sys.stderr, "usage: " + sys.argv[0] + " -d <gff-folder>"
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -p        pfam annotation of proteins" )
  stdout( " -m        misa file incl. protein in last column" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hm:p:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-m': args['misa'] = value
    if key == '-p': args['pfam'] = value
    
  if not args.has_key('misa'):
    print >> sys.stderr, "misa file argument missing."
    show_help()
  elif not file_exists( args.get('misa') ):
    print >> sys.stderr, "misa file does not exist."
    show_help()

  if not args.has_key('pfam'):
    print >> sys.stderr, "pfam file argument missing."
    show_help()
  elif not file_exists( args.get('pfam') ):
    print >> sys.stderr, "pfam file does not exist."
    show_help()

  return args


# =============================================================================
def get_ssrs(file):
    hash = defaultdict(list)
    fo = open(file)
    for line in fo: 
      if line.startswith("ID\t"): continue
      m = MisaSSR(line)
      hash[m.geneid].append(m)
    fo.close()
    print >> sys.stderr, "read %s microsatellites" % len(hash)
    return hash

# =============================================================================
def get_pfam(file):
  hash = defaultdict(list)
  counthash = defaultdict(int)
  fo = open(file)
  for line in fo:
    pd = PfamDomain(line)
    hash[pd.get_attr('seq_id')].append(pd)
    counthash[pd.get_attr('hmm_acc') + '|' + pd.get_attr('hmm_name')] += 1
  fo.close()
  print >> sys.stderr, "read %s pfam annotations" % len(hash)
  return hash, counthash

# =============================================================================
def sd(pylist):
  rvec = robjects.FloatVector(pylist)
  sd = R['sd'](rvec)[0]
  return sd

# =============================================================================
def pnorm(value, mean, sd):
  p = R['pnorm'](value, mean=mean, sd=sd, lower=False)[0]
  return p

# =============================================================================
def fisher_test(pylist):
  rcountsvec = robjects.IntVector(pylist)
  rmatrix = R['matrix'](rcountsvec,2,2)
  p = R['fisher.test'](rmatrix, alternative="greater")[0][0]
  p = float(p)
  pylist.append(p)
  print >> sys.stderr, string.join([str(x) for x in pylist],"\t")
  return p

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  SSRs = get_ssrs(args['misa']) 
  PFAMs, pfamglobalcounts = get_pfam(args['pfam']) 
  pfamhits = defaultdict(int)
  for sid, SSRs in SSRs.iteritems():
    for SSR in SSRs:
      for PFAM in PFAMs[SSR.feature]:
        if SSR.startpos < int(PFAM.get_attr('alignment_start')) and SSR.endpos < int(PFAM.get_attr('alignment_start')): continue
        if SSR.startpos > int(PFAM.get_attr('alignment_end')) and SSR.endpos > int(PFAM.get_attr('alignment_end')): continue
        pfamhits[PFAM.get_attr('hmm_acc') + '|' + PFAM.get_attr('hmm_name')] += 1

  totalssrsindomains = 0
  for pfam, count in pfamhits.iteritems(): totalssrsindomains += count
  totaldomains = 0
  for id, count in pfamglobalcounts.iteritems(): totaldomains += count

  ratios = {}
  totalssrcount = sum(pfamhits.values())
  totaldomaincount = sum(pfamglobalcounts.values())
  expectedfreq = 1.0*totalssrcount/totaldomaincount
  for id, count in pfamhits.iteritems():
    p = fisher_test([count, max([0,pfamglobalcounts[id] - count]), totalssrcount - count, totaldomaincount - max([0,pfamglobalcounts[id] - count])])
    ratios[id] = p

  pvalues = []
  ids = []
  for id, p in ratios.iteritems(): 
    pvalues.append(p)
    ids.append(id)

  ps = robjects.FloatVector(pvalues)
  psadjusted = tuple(R['p.adjust'](ps, method="fdr"))
  for i in range(len(psadjusted)):
    if psadjusted[i] < 0.05:
      print ids[i] + "\t" + str(ratios[ids[i]]) + "\t" + str(pvalues[i]) + "\t" + str(psadjusted[i])

#  for id, p in testeddomains.iteritems():
#    if p < 0.05: print id + "\t" + str(p)
#    print id + "\t" + str(p) #+ "\t" + str(pfamglobalcounts[id] - count) + "\t" + str(totalssrsindomains - count) + "\t" + str(totaldomains - totalssrsindomains - pfamglobalcounts[id] - count)



# =============================================================================
args = handle_arguments()
main( args )
