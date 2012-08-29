import sys

class OrthoMCLCluster():
  def __init__(self, line):
    descr, genedefs = line.split("\t")
    genedefs = genedefs.split()
    self.name = descr[:descr.index('(')].lower()
    self.geneHash = {}
    self.speciesHash = {}
    for genedef in genedefs:
      geneid = genedef[:genedef.index('(')]
      species = genedef[genedef.index('(')+1:-1]
      self.geneHash[geneid] = species
      if self.speciesHash.has_key(species): self.speciesHash[species].append(geneid)
      else: self.speciesHash[species] = [geneid]

  def add_gene(self, geneid, species):
    if not self.geneHash.has_key(geneid):
      self.speciesHash[species].append(geneid)
      self.geneHash[geneid] = species
  def get_name(self): return self.name
  def get_count(self): return len(self.geneHash)
  def get_gene_hash(self): return self.geneHash
  def get_species_hash(self): return self.speciesHash
  def to_s(self):
    sys.stdout.write(self.name + "(" + str(len(self.geneHash)) + " genes, " + str(len(self.speciesHash)) + ")\t")
    first = 1
    for geneid, species in self.geneHash.iteritems():
      if first == 0: sys.stdout.write(" ")
      first = 0
      sys.stdout.write(geneid + "(" + species + ")")
    sys.stdout.write("\n")

