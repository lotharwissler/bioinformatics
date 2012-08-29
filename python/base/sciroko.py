import string

class SSR():
  def __init__(self, line):
    columns = line.rstrip().split("\t")
    self.seqname = columns.pop(0)
    self.motif = columns.pop(0)
    self.motif_std = columns.pop(0)
    self.startpos = int(columns.pop(0))
    self.endpos = int(columns.pop(0))
    self.length = int(columns.pop(0))
    self.score = int(columns.pop(0))
    self.mismatches = int(columns.pop(0))
    self.mode = columns.pop(0)
    self.seq = columns.pop(0)

  def to_s(self):
    array = [self.seqname, self.motif, self.motif_std, str(self.startpos), str(self.endpos), str(self.length), str(self.score), str(self.mismatches), self.mode, self.seq]
    return string.join(array, "\t")

  def is_perfect(self):
    if self.mismatches == 0: return 1
    return 0

  def motif_len(self): return len(self.motif_std)

  def compare_to(self, other):
    if self.organism == other.organism or self.chromosome != other.chromosome: return 0 # not comparable
    if self.motif_std != other.motif_std: return 99 # different motif, i.e. no similarity between the SSRs
    if self.is_perfect and other.is_perfect:
      if self.seq == other.seq: return 1 # 100% perfect matches
      if self.length != other.length: return 11 # two perfect SSRs with polymorphic change

    if (self.is_perfect and not other.is_perfect) or (not self.is_perfect and other.is_perfect):
      if self.length == other.length: return 31

    if not self.is_perfect and not other.is_perfect:
      if self.seq == other.seq: return 51


  def is_perfect_match_to(self, other):
    if self.motif_std != other.motif_std: return 0
    if self.length != other.length: return 0
    if self.mismatches != other.mismatches: return 0
    return 1

  def is_imperfect_match_to(self, other):
    if self.motif_std != other.motif_std: return 0
    if self.length != other.length: return 0
    return 1

  def is_polymorphic_to(self, other):
    if self.motif != other.motif: return 0
    if self.repeats == other.repeats: return 0
    return 1

  def is_shifted_to(self, other):
    if self.motif == other.motif: return 0
    if self.type != other.type: return 0
    m = self.motif
    for i in range(len(self.motif)):
      m = m[1:] + m[0]
      if m == other.motif: return 1
    return 0


class MisaSSR():
  def __init__(self, line):
    self.feature = 0
    columns = line.rstrip().split("\t")
    self.geneid = columns.pop(0)
    self.ssrnr = int(columns.pop(0))
    self.type = columns.pop(0)
    self.pattern = columns.pop(0)
    self.length = int(columns.pop(0))
    self.startpos = int(columns.pop(0))
    self.endpos = int(columns.pop(0))
    if len(columns) > 0: self.feature = columns.pop(0)
    if self.type != "c" and self.type != "c*":
      self.motif = self.pattern[1:self.pattern.index(")")]
      if self.pattern.endswith("*"): self.repeats = int(self.pattern[self.pattern.index(")")+1:-1])
      else: self.repeats = int(self.pattern[self.pattern.index(")")+1:])

  def to_s(self):
    array = [self.geneid, str(self.ssrnr), self.type, self.pattern, str(self.length), str(self.startpos), str(self.endpos)]
    return string.join(array, "\t")

  def is_perfect_match_to(self, other):
    if self.pattern != other.pattern: return 0
    return 1

  def is_polymorphic_to(self, other):
    if self.motif != other.motif: return 0
    if self.repeats == other.repeats: return 0
    return 1

  def is_shifted_to(self, other):
    if self.motif == other.motif: return 0
    if self.type != other.type: return 0
    m = self.motif
    for i in range(len(self.motif)):
      m = m[1:] + m[0]
      if m == other.motif: return 1
    return 0


