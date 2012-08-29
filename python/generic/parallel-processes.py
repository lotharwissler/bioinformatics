#!/usr/bin/python
import sys, os
import math, time
import threading
import getopt         # comand line argument handling
from low import *

# =============================================================================
def show_help():
  """ displays the program parameter list and usage information """
  stdout( "usage: " + sys.argv[0] + " -f <path> -n " )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        file that contains all system calls, one per line" )
  stdout( " -n        number of processes to be run in parallel" )
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

  args = {}
  for key, value in keys:
    if key == '-f': args['file'] =  value
    if key == '-n': args['ncpu'] =  int(value)
    
  if not args.has_key('file'):
    stderr( "process file missing." )
    show_help()
  if not file_exists( args.get('file') ):
    stderr( "process file does not exist." )
    show_help()

  if not args.has_key('ncpu'):
    stderr( "number of CPUs to use is missing." )
    show_help()

  return args



# =============================================================================
def get_jobs( file ):
  jobs = []
  fo = open( file )
  for line in fo:
    if not line.startswith("#"): jobs.append(line.rstrip())
  fo.close()
  return jobs

# =============================================================================
class MyThread( threading.Thread ):
  
  def set_command(self, command):
    self.command = command

  def run(self):
    os.system(self.command)


# =============================================================================
# =============================================================================
def main( args ):
  
  jobs = get_jobs( args.get('file') )
  totaljobs = len(jobs)
  infomsg( "Collected %s jobs queued | will distribute them among %s CPUs" %(len(jobs), args.get('ncpu')) )
  
  start_time = time.time()
  while threading.activeCount() > 1 or len(jobs) > 0:
    # check existing threads: still running?
    # fill up all remaining slots
    elapsed = time.time() - start_time
    while threading.activeCount() <= args.get('ncpu') and len(jobs) > 0:
      # start new thread
      cmd = jobs.pop(0)
      t = MyThread()
      t.set_command( cmd )
      t.start()

    remain = elapsed / (totaljobs - len(jobs) + (threading.activeCount() -1)) * len(jobs)
    info( "\telapsed: %s\tremaining: %s\t[ jobs done: %s | remain: %s | active: %s ] " % (humanize_time(elapsed), humanize_time(remain), totaljobs - len(jobs) - (threading.activeCount() -1), len(jobs), threading.activeCount() -1) )
    time.sleep(0.2)
    
  info( "\n" )
  infomsg( "DONE." )


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
