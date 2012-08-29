#!/usr/bin/python

import sys,os,getopt

OUTFILEPART2 = 'tmp.orf.part2.fasta'

#==============================================================================
def usage():
  print """Hello!
  
  Following options are possible:
  -i:\tparsed BLASTX best hit definitions
  -j:\tinput sequences in FASTA format
  -t:\tminimum length for in silico predicted ORFs
  """

#==============================================================================
def main( XMLfile, CAP3, threshold ):
  # First, elongate BLAST-hits

  os.system("orf_prediction_part1.py -b "+str(XMLfile)+" -f "+str(CAP3) )
  #print "BLASTelongator has finished. Starting 2nd part..."
  # It has written to temp and now comes Ina's script

  os.system("orf_prediction_part2.py"+" -t "+str(threshold)+" -f " + OUTFILEPART2 )
  #print "ORF-Prediction has finished. Removing temp-files.."
  #os.system("cat BLASTelongatorHits.out SimulatedORFS.out > "+str(outfile))  
  #os.system("rm BLASTelongatorHits.out")
  os.system("rm " + OUTFILEPART2)
  #os.system("rm SimulatedORFS.out")
  #print "Done. See you soon!"



#==============================================================================
# MAIN ========================================================================
#==============================================================================
try:  
  if len(sys.argv) > 1:
    opts, args = getopt.getopt(sys.argv[1:],"i:j:t:h")
  else:
    usage()
    #print "Hello! I take at least 3 arguments. I have the following options: -i defines the input XML-file which you want to use, -j defines the CAP3-outputfile, -t defines the threshold for in silico predicted proteins, -o defines the outfile, -h gives you more help! See you soon! You provided:" 
    sys.exit()
except getopt.GetoptError, err:
  print "Something went wrong - maybe this helps: " + str(err)
  sys.exit()

for o, a in opts:
  if o == "-h":
    usage()
    sys.exit()
  elif o == "-i":
    if os.path.exists(a):
      if os.path.isfile(a):
        XMLfile = a
      elif os.path.isdir(a):
        print "Specified XML-file is a directory!"
    else:
      print "Something is wrong with the XML-file, maybe it doesn't exist?"
  elif o == "-j":
    if os.path.exists(a):
      if os.path.isfile(a):
        CAP3 = a
      elif os.path.isdir(a):
        print "Specified CAP3-file is a directory!"
    else:
      print "Something is wrong with the CAP3-file, maybe it doesn't exist?"
  elif o == "-t":
    threshold = a
  else:
    print "Something went wrong ;_;. Maybe the file you specified doesn't exist?"

if len(opts) == 3:
  main( XMLfile, CAP3, threshold )
else:
  print len(opts)
  print "Again, hello to you! You do not have the required amount of arguments given. Please specify them. For more, see -h! I AM THE PREDICTOR"
