#!/usr/bin/python

import os, sys     # low level handling, such as command line stuff
import string      # string methods available
import re          # regular expressions
import getopt      # comand line argument handling
from low import *  # custom functions, written by myself


# =============================================================================  
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        signalp output file (short format)" )
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

  args = {}
  for key, value in keys:
    if key == '-f': args['signalpfile'] = value
    
  for key in ['signalpfile']:
    if key.endswith("file"):
      if not args_file_exists(args, key): show_help()
    elif key.endswith("dir"):
      if not args_dir_exists(args, key): show_help()
    elif not args.has_key(key):
      print >> sys.stderr, "missing argument", key
      show_help()
  return args

# =============================================================================
def statusbar(current, total, message="", width=40):
  progress = 1.0*current/total
  if message != "": message = "[" + message + "]"
  progressbar = "=" * int(progress*width)
  while len(progressbar) < width: progressbar += " " 
  sys.stderr.write("\r   0% " + progressbar + " 100% " + message)
  if progress == 1.0: sys.stderr.write("\n")
  
# =============================================================================
# name                Cmax  pos ?  Ymax  pos ?  Smax  pos ?  Smean ?  D     ? 	# name      !  Cmax  pos ?  Sprob ?
class SignalpResult:
  def __init__(self, line):
    col = line.rstrip().split()
    self.gid = col.pop(0)
    self.NN_Cmax = float(col.pop(0))
    self.NN_Cmax_pos = int(col.pop(0))
    if col.pop(0) == 'Y': self.NN_Cmax_sign = True
    else: self.NN_Cmax_sign = False
    self.NN_Ymax = float(col.pop(0))
    self.NN_Ymax_pos = int(col.pop(0))
    if col.pop(0) == 'Y': self.NN_Ymax_sign = True
    else: self.NN_Ymax_sign = False
    self.NN_Smax = float(col.pop(0))
    self.NN_Smax_pos = int(col.pop(0))
    if col.pop(0) == 'Y': self.NN_Smax_sign = True
    else: self.NN_Smax_sign = False
    
    self.NN_Smean = float(col.pop(0))
    if col.pop(0) == 'Y': self.NN_Smean_sign = True
    else: self.NN_Smean_sign = False
    self.NN_D = float(col.pop(0))
    if col.pop(0) == 'Y': self.NN_D_sign = True
    else: self.NN_D_sign = False
    
    col.pop(0)
    self.HMM_res = col.pop(0)
    self.HMM_Cmax = float(col.pop(0))
    self.HMM_Cmax_pos = int(col.pop(0))
    if col.pop(0) == 'Y': self.HMM_Cmax_sign = True
    else: self.HMM_Cmax_sign = False
    self.HMM_Sprob = float(col.pop(0))
    if col.pop(0) == 'Y': self.HMM_Sprob_sign = True
    else: self.HMM_Sprob_sign = False
    
  def is_significant(self, NN=True, HMM=True):
    sign = True
    if not NN and not HMM:
      return [self.NN_Cmax_sign, self.NN_Ymax_sign, self.NN_Smax_sign, self.NN_Smean_sign, self.NN_D_sign, self.HMM_Cmax_sign, HMM_Sprob_sign].count(True)
    else:
      if NN:
        if not self.NN_Cmax_sign and not self.NN_Ymax_sign and not self.NN_Smax_sign and not self.NN_Smean_sign and not self.NN_D_sign: sign = False
      if HMM: 
        if not self.HMM_Cmax_sign and not self.HMM_Sprob_sign: sign = False
      return sign
    
    
    
# =============================================================================
# === MAIN ====================================================================
# =============================================================================
def main( args ):
  fo = open(args["signalpfile"])
  for line in fo:
    if line.startswith("#"): continue
    #print line
    sr = SignalpResult(line)
    if sr.is_significant(True, True): print sr.gid
  fo.close()
  

# =============================================================================
args = handle_arguments()
main( args )

