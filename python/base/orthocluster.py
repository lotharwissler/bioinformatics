import sys

# =============================================================================
def parse(clusterfile, poshash={}):
  nspecies = 0
  clusterhash = {}
  fo = open(clusterfile)
  while 1:
    line = fo.readline()
    if not line: break
    if line.startswith("No. of sequence"): nspecies = int(line.split()[-1])
    if line.startswith("CL-"):
      cols = line.split()
      sc = OrthoCluster(cols[0])
      ngenes = int(max(cols[1:nspecies+1]))
      fo.readline()
      for i in range(ngenes):
        cols = fo.readline().split()
        for j in range(nspecies):
          cols.pop(0)
          cols.pop(0)
          strand = cols.pop(0)
          scaffold = cols.pop(0)
          geneid = cols.pop(0)
          if poshash.has_key(geneid):
            startpos, endpos = poshash[geneid][1:3]
          else:
            print >> sys.stderr, "geneid", geneid, "not found in poshash"
            startpos, endpos = None, None
          sr = SyntenicRegion(geneid, scaffold, strand, startpos, endpos)
          sc.add_syntenic_region(sr, j)
      clusterhash[sc.id] = sc
  fo.close()
  return nspecies, clusterhash


# =============================================================================
class SyntenicRegion():
  def __init__(self, geneid, scaffold, strand, startpos, endpos):
    self.geneid = geneid
    self.scaffold = scaffold
    self.strand = strand
    self.startpos = startpos
    self.endpos = endpos

# =============================================================================
class OrthoCluster():
  def __init__(self, clusterid):
    self.id = clusterid
    self.syntenic_regions = {}

  def add_syntenic_region(self, sr, index):
    if not self.syntenic_regions.has_key(index): self.syntenic_regions[index] = []
    self.syntenic_regions[index].append(sr)
