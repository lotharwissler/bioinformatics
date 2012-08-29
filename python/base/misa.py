import string

class MisaSSRspecies():
  def __init__(self, line):
    self.feature = 0
    columns = line.rstrip().split("\t")
    self.species = columns.pop(0)
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
    array = [self.species, self.geneid, str(self.ssrnr), self.type, self.pattern, str(self.length), str(self.startpos), str(self.endpos)]
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


