import sys, string

def get_gff_hash(gffile):
  hash = {}
  fo = open(gffile)
  for line in fo:
    gf = GeneFeature(line)
    if not hash.has_key(gf.seqid): hash[gf.seqid] = []
    hash[gf.seqid].append(gf)
  fo.close()
  return hash
  

class GeneFeature():
  def __init__(self, line):
    columns = line.rstrip().split("\t")
    if not len(columns) == 9:
      print >> sys.stderr, "GFF3 with incorrect number of columns. Expected: 9 | Observed: %s" % len(columns)
      print >> sys.stderr, "\"%s\"" % line
      sys.exit(1)
    self.seqid = columns.pop(0)
    self.source = columns.pop(0)
    self.ftype = columns.pop(0)
    self.start = int(columns.pop(0))
    self.stop = int(columns.pop(0))
    self.score = columns.pop(0)
    self.strand = columns.pop(0)
    self.phase = columns.pop(0)
    self.attributes = columns.pop(0)

  def get_attributes(self):
    hash = {}
    for e in self.attributes.split(";"):
      if e == '': continue
      k, v = e.split("=")
      hash[k] = v
    return hash
    
  def set_attribute(self, key, value):
    hash = {}
    for e in self.attributes.split(";"):
      if e == '': continue
      k, v = e.split("=")
      hash[k] = v
    if hash.has_key(key):
      hash[key] = value
      self.attributes = ""
      for k, v in hash.iteritems(): self.attributes += "%s=%s;" %(k, v)
    else:
      self.attributes += "%s=%s;" %(key, value)
    

  def to_string(self):
    return string.join([self.seqid, self.source, self.ftype, str(self.start), str(self.stop), self.score, self.strand, self.phase, self.attributes], "\t")
