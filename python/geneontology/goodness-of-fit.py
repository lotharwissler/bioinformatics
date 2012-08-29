#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
import sqlite3
import glob
import newick
import pylab
import rpy2.robjects as rpy2
import numpy
import copy
from low import *  # custom functions, written by myself


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        CSV file with two follows (first x, second y)" )
  stdout( " -t        type of fit [linear|log|exp]" )
  stdout( " " )
  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the data dir """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()  
  
  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:t:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = {'fit': "linear"}
  for key, value in keys:
    if key == '-f': args['csvfile'] = value
    if key == '-t': args['fit'] = value
    
  if not args.has_key('csvfile'):
    stderr( "csv file argument missing." )
    show_help()
  elif not file_exists( args.get('csvfile') ):
    stderr( "csv file does not exist." )
    show_help()
  
  return args


# =============================================================================
def statusbar(counter, total, message="", width=50):
  fraction = 1.0*counter/total
  progressbar = int(width*fraction) * "="
  while len(progressbar) < width: progressbar += " "
  sys.stderr.write("\r 0% [" + progressbar + "] 100% ")
  if message != "": sys.stderr.write("| " + message)
  if fraction == 1.0: sys.stderr.write("\n")

# =============================================================================
def get_xy(csvfile):
  x, y = [], []
  fo = open(csvfile)
  for line in fo:
    col = line.rstrip().split(',')
    x.append(float(col[0]))
    y.append(float(col[1]))    
  fo.close()
  return x, y
  

# =============================================================================
def get_fit_params(x, y, fit):
  R = rpy2.r
  rpy2.globalEnv['x'] = rpy2.FloatVector(x)
  rpy2.globalEnv['y'] =  rpy2.FloatVector(y)
  ymean = numpy.average(y)
  SStot, SSerr = [], []
  
  if fit == 'linear':
    lmfit = R.lm('y ~ x')
  elif fit == 'log':
    lmfit = R.lm('y ~ log(x)')
  elif fit == 'exp':
    lmfit = R.lm('y ~ exp(x)')
  a, b = lmfit[0][0], lmfit[0][1]
  
  for i in range(len(x)):
    SStot.append(numpy.square(y[i]-ymean))
    if fit == 'linear':
      SSerr.append(numpy.square(y[i]-(a+x[i]*b)))
    elif fit == 'log':
      SSerr.append(numpy.square(y[i]-(a+numpy.log(x[i])*b)))
    elif fit == 'exp':
      SSerr.append(numpy.square(y[i]-(a+numpy.exp(x[i])*b)))
  SStot, SSerr = sum(SStot), sum(SSerr)
  rsquared = 1 - (SSerr / SStot)
  return a, b, rsquared
  

# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  x, y = get_xy(args['csvfile'])
  a, b, rsquared = get_fit_params(x, y, args['fit'])
  print rsquared


# =============================================================================
args = handle_arguments()
main( args )


