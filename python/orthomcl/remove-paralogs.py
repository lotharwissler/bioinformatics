#!/usr/bin/python
import os, sys
import warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)

# takes an input protein fasta file and an orthomcl.gg file
# orthomcl.gg file format:
# speciesname1: id1 id2 id3 id4 .... full genome
# speciesname2: id1 id2 id3 id4 .... full genome
#
# with these infos, the goal is to get only one protein sequence per species
# we use t-coffee to find the most similar protein sequence per species
# to the whole cluster. so in case one species contributes several sequences 
# to a cluster, we choose the one species to keep which has the highest average 
# similarity to the rest of the cluster. if more than 1 sequence yield the highest
# avgsim, we determine whether these protein sequences are (1) all identical, 
# or whether they are (2) slightly different. In case (1), we choose any sequence
# randomly because it does not matter. In case (2), we sum up all pairwise
# similarities for each candidate sequence, and keep only the one sequence
# with the highest sum. If these are identical as well, we again choose randomly
# (should happen very rarely).



def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " fasta-folder  all-proteins.fasta  orthomcl.gg"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 4: usage()
  inFolder, allProteins, inGG = sys.argv[1:5]
  if not inFolder.endswith("/"): inFolder += '/'
  return inFolder, allProteins, inGG


def get_number_of_species(inGG):
  count = 0
  fo = open(inGG)
  for line in fo: count += 1
  fo.close()
  return count


def read_gg(inGG):
  outHash = {}
  speciesArray = []
  fo = open(inGG)
  for line in fo:
    line = line.rstrip()
    cols = line.split()
    species = str(cols[0])[:-1]
    if not species in speciesArray: speciesArray.append(species)
    for col in cols[1:]:
      outHash[col] = species
  fo.close()
  return outHash, speciesArray


def parse_sim_out(inSeqs,inSim):
  allIDs = []
  fo = open(inSeqs)
  for line in fo:
    if line.startswith(">"): 
      line = line.rstrip()
      allIDs.append(line[1:])
  fo.close()

  def translate_id(id, aI=allIDs):
    if not id in aI:
      for E in aI: 
        e = E
        e = e.replace("#","_")
        e = e.replace(":","_")
        if id == e: id = E
    return id

  outHash = {}
  fo = open( inSim )
  for line in fo:
    line = line.rstrip()
    cols = line.split()
    if line.startswith("AVG"):
      id, avgsim = translate_id(cols[2]), float(cols[4])
      outHash[id].insert(0, avgsim)
    elif line.startswith("TOP") or line.startswith("BOT"):
      id1, id2, sim = translate_id(cols[4]), translate_id(cols[5]), float(cols[6])
      outHash[id1 + "$$$" + id2] = sim
      if outHash.has_key(id1): outHash[id1].append(sim)
      else: outHash[id1] = [sim]
  return outHash


def tcoffee_similarity(inSeqs, inSpeciesHash):
  out = inSeqs + ".tcoffee.sim"
  if not os.path.exists(out) or not os.path.getsize(out) > 0:
    syscall = '~/bin/t-coffee -other_pg seq_reformat -in %s -output sim 1> %s' %(inSeqs, out)
    ec = os.system(syscall)
    if not ec == 0 or not os.path.isfile(out) or os.path.getsize(out) == 0:
      print "T-COFFEE did not run smoothly. check manually: " + syscall
      sys.exit(3)
  
  simHash = parse_sim_out(inSeqs,out)
  specHash = {}
  for id, sim in simHash.iteritems():
    if id.count('$$$') > 0: continue
    if not inSpeciesHash.has_key(id):
      id = id.replace("___","_#_")
    species = inSpeciesHash[id]
    if not specHash.has_key(species) or specHash[species][0] < sim[0]:
      specHash[species] = [sim[0], id]
    elif specHash.has_key(species) and specHash[species][0] == sim[0]: specHash[species].append(id)


  idsToKeep = []
  for species, array in specHash.iteritems():
    avgsim = array[0]
    ids = array[1:]
    while len(ids) > 1:
      id1, id2 = ids[0], ids[1]
      sim = simHash[id1 + '$$$' + id2]
      if sim == 100.00: ids.pop(1)
      else: 
        simid1 = sum(simHash[id1])
        simid2 = sum(simHash[id2])
        if simid1 == simid2:
          print "WARNING %s: identical avgsims but no identity between sequences %s [%s] / %s [%s] (%s | %s)" %(out, id1, simid1, id2, simid2, sim, avgsim)
          ids.pop(1)
        elif simid1 > simid2: ids.pop(1)
        elif simid1 < simid2: ids.pop(0)
    idsToKeep.append(ids[0])
      
  return idsToKeep


def cache_genomes(file, recreate=0):
  outdbm = file + ".dbm"
  if os.path.exists(outdbm) and os.path.getsize(outdbm) > 0 and not recreate: return outdbm
  DBM = anydbm.open( outdbm, 'c' )
  fo = open(file)
  key = ""
  for line in fo: 
    line = line.strip()
    if line.startswith(">"): 
      key = line[1:]
      if key.count(" ") > 0: key = key[:key.index(" ")]
      DBM[key] = ""
    else:
      DBM[key] += line
  DBM.close()
  return outdbm

  

def slim_cluster(inProteins, seqHash, speciesHash, speciesArray):
  outdir = "paralog-free-clusters/"
  if not os.path.exists(outdir): os.mkdir(outdir)
  sout, serr = catch_bash_cmd_output( "grep '>' -c %s" %(inProteins) )
  nGenes = int( sout )
  outbase = os.path.split(inProteins)[1]
  outFasta = outdir + os.path.splitext(outbase)[0] + ".fasta"
  outUfasta = outdir + os.path.splitext(outbase)[0] + ".ufasta"
  
  if nGenes > len(speciesArray):
    idsToKeep = tcoffee_similarity(inProteins, speciesHash)
  else:
    idsToKeep = []
    fo = open(inProteins)
    for line in fo:
      if line.startswith(">"): 
        line = line.rstrip()
        idsToKeep.append(line[1:])
    fo.close()

  outstring = os.path.splitext(outbase)[0] + "(" + str(len(speciesArray)) + " genes, " + str(len(speciesArray)) + " taxa):\t"
  fwi = open(outFasta, 'w')
  fwu = open(outUfasta, 'w')
  for id in idsToKeep: 
    fwi.write(">" + id + "\n")
    fwu.write(">" + speciesHash[id] + "\n")
    sequence = seqHash[id]
    i = 0
    width = 60
    while i < len(sequence):
      frac = sequence[i:min([len(sequence),i+width])]
      fwi.write(frac + "\n")
      fwu.write(frac + "\n")
      i += width
    outstring += id + "(" + speciesHash[id] + ") "
  fwi.close()
  fwu.close()
  return outstring
 

def get_all_fasta_files(inFolder):
  files = []
  for inFile in os.listdir(inFolder):
    if not inFile.endswith(".fasta"): continue
    files.append(inFolder + inFile)
  return files


def main():
  inFolder, allProteins, inGG = plausi()
  inFiles = get_all_fasta_files(inFolder)
  speciesHash, speciesArray = read_gg(inGG)
  dbm = cache_genomes(allProteins)
  seqHash = anydbm.open(dbm, "r")
  total, count = len(inFiles), 0
  fw = open("noparalogs.orthomcl.out", "w")
  for inFile in inFiles:
    outstring = slim_cluster(inFile, seqHash, speciesHash, speciesArray)
    fw.write( outstring + "\n" )
    count += 1
    progress = int(50.0 * count / total) * "#"
    progress += (50 - len(progress)) * " "
    info("       0% " + progress + " 100%     ")
  info("       0% " + progress + " 100%    \n")
  seqHash.close()
  fw.close()


main()
