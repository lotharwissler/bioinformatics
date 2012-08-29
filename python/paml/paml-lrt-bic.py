#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
from low import *			# custom functions, written by myself
import math

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        nt alignment file" )
  stdout( " -s        simple mode: only test between M0 and FreeRatio" )
  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:t:s" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {}
  for key, value in keys:
    if key == '-f':	args['file'] = value
    if key == '-s':	args['simple'] = 1
        
  if not args.has_key('file'):
    stderr( "file missing." )
    show_help()
  if not file_exists( args.get('file') ):
    stderr( "file does not exist." )
    show_help()
    
  return args


# =============================================================================
def LRT( lnL1, lnL2 ):
  """ 
  The LRT statistic, or twice the log likelihood difference between the two compared models (2{Delta}{ell}), may be compared against {chi}Formula, with critical values to be 5.99 and 9.21 at 5% and 1% significance levels, respectively.
  """
  lnL1 = float(lnL1)
  lnL2 = float(lnL2)
  return (2 * math.fabs( (lnL1 - lnL2) ))

# =============================================================================
def BIC( np, lnL, length ):
  """
  BIC = -2 ln( L ) + k ln(n),
  n = number of observations = sample size = alignment length (nt)
  k = number of parameters
  ln( L ) = PAML lnL = log likelihood
  """
  lnL = float(lnL)
  np = int(np)
  length = int(length)
  return ( ((-2)* lnL) + (np * math.log(length)) )

# =============================================================================
# =============================================================================
def main( args ):
  
  fo = open( args.get('file'), 'r' )
  for line in fo:
    line = line.rstrip()
    columns = line.split("\t")
    name = columns.pop(0)
    length = columns.pop(0)
    MH = {}
    while len( columns ) > 2:
      hash = {}
      Modelname = columns.pop(0)
      np = int(columns.pop(0))
      lnL = float(columns.pop(0))
      MH[ Modelname ] = [ np, lnL, Modelname ]
    
    if args.has_key('simple'):
      # Free vs. M0
      lrt = LRT( MH.get("Free")[1], MH.get("M0")[1])
      df = 2* ( int(MH.get("Free")[0]) - int(MH.get("M0")[0]) )
      L = [name, str(lrt), str(df)]
      print string.join( L, "\t" )

    else:
      # M7 vs. M8
      lnL_M7 = MH.get("M7")[1]
      lnL_M8 = MH.get("M8")[1]
      if lnL_M7 > lnL_M8: MH[ "M7M8" ] = MH.get("M7")
      else: MH[ "M7M8" ] = MH.get("M8")

      # M3K2 vs. M0
      if MH.get("M3K2")[1] > MH.get("M0")[1]: MH[ "M3K2M0" ] = MH.get("M3K2")
      else: MH[ "M3K2M0" ] = MH.get("M0")
      
      # M3K3 vs. M0
      if MH.get("M3K3")[1] > MH.get("M0")[1]: MH[ "M3K3M0" ] = MH.get("M3K3")
      else: MH[ "M3K3M0" ] = MH.get("M0")

      # M3 vs. M0
      if MH.get("M3K2")[1] > MH.get("M0")[1] and MH.get("M3K3")[1] > MH.get("M0")[1]:
        # BIC
        BIC_M3K2 = BIC( MH.get("M3K2")[0],  MH.get("M3K2")[1], length )
        BIC_M3K3 = BIC( MH.get("M3K3")[0],  MH.get("M3K3")[1], length )
        if BIC_M3K2 < BIC_M3K3: MH[ "M3M0" ] = MH.get("M3K2")
        else: MH[ "M3M0" ] = MH.get("M3K3")
      elif MH.get("M3K2")[1] > MH.get("M0")[1]: MH[ "M3M0" ] = MH.get("M3K2")
      elif MH.get("M3K3")[1] > MH.get("M0")[1]: MH[ "M3M0" ] = MH.get("M3K3")
      else: MH[ "M3M0" ] = MH.get("M0")

      # Free vs. M0
      if MH.get("Free")[1] > MH.get("M0")[1]: MH[ "FreeM0" ] = MH.get("Free")
      else: MH[ "FreeM0" ] = MH.get("M0")

      # now compare winners of
      # - M0 vs Free
      # - M0 vs M3.2/3
      # - M7 vs M8

      # M3M0 vs. FreeM0
      HB = {}
      for M in ["M3M0", "M7M8", "FreeM0"]:
        HB[ MH.get(M)[2] ] = BIC( MH.get(M)[0], MH.get(M)[1], length )
      
      for name, array in MH.iteritems():
        print name, "--", array
#      print "MH:", MH
      print "HB:", HB
    
      sortedKeys = sort_by_value( HB )
      L = []
      for key in sortedKeys[0:1]:
        L.append( key )
        L.append( str(HB.get(key)) )
      print string.join( L, "\t" )
      #BIC_M3M0 = BIC( MH.get("M3M0")[0],  MH.get("M3M0")[1], length )
      #BIC_FreeM0 = BIC( MH.get("FreeM0")[0],  MH.get("FreeM0")[1], length )
      #BIC_M7M8 = BIC( MH.get("M7M8")[0],  MH.get("M7M8")[1], length )
    
    

# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
