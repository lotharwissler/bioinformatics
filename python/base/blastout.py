import string

# =============================================================================
class BlastHit:
  def __init__(self, line):
    cols = line.split("\t")
    self.qid, self.hid = cols.pop(0), cols.pop(0)
    self.identity = float(cols.pop(0))
    self.alnlen = int(cols.pop(0))
    self.mismatch = int(cols.pop(0))
    self.gap = int(cols.pop(0))
    self.qstart = int(cols.pop(0))
    self.qstop = int(cols.pop(0))
    self.hstart = int(cols.pop(0))
    self.hstop = int(cols.pop(0))
    self.evalue = float(cols.pop(0))
    self.score = float(cols.pop(0))
    
  def to_s(self):
    out = []
    out += [self.qid, self.hid, str(self.identity), str(self.alnlen)]
    out += [str(self.mismatch), str(self.gap), str(self.qstart), str(self.qstop)]
    out += [str(self.hstart), str(self.hstop), str(self.evalue), str(self.score)]
    return string.join(out, "\t")

# =============================================================================
def get_query_hash(blastoutfile, evalue=10.0):
  qh = {}
  fo = open(blastoutfile)
  for line in fo:
    line = line.rstrip()
    if len(line) == 0 or line.startswith('#') or not len(line.split("\t")) == 12: continue
    blasthit = BlastHit(line)
    if blasthit.evalue > evalue: continue
    if not qh.has_key(blasthit.qid): qh[blasthit.qid] = []
    qh[blasthit.qid].append(blasthit)
  fo.close()
  return qh

# =============================================================================
def get_sequence_hash(fastafile):
  seqhash = {}
  key = ""
  fo = open(fastafile)
  for line in fo:
    if line.startswith(">"):
      gid = re.match(">(\S+)", line).group(1)
      key = gid
      seqhash[key] = ""
    else:
      if key != "": seqhash[key] += line.strip()
  fo.close()
  return seqhash
