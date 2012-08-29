class GOTerm():
  def __init__(self, lines):
    self.id = ""
    self.name = ""
    self.namespace = ""
    self.definition = ""
    self.is_a = []
    self.alt_ids = []
    self.xrefs = []
    self.synonyms = []
    self.obsolete = 0

    for line in lines:
      line = line.strip()
      if line.startswith("id: "): self.id = line[line.index(":")+2:]
      if line.startswith("name: "): self.name = line[line.index(":")+2:]
      if line.startswith("namespace: "): self.namespace = line[line.index(":")+2:]
      if line.startswith("def: "): self.definition = line[line.index(":")+2:]
      if line.startswith("is_a: "): self.is_a.append( line[line.index(":")+2:] )
      if line.startswith("is_obsolete: true"): self.obsolete = 1
      if line.startswith("alt_id: "): self.alt_ids.append( line[line.index(":")+2:] )
      if line.startswith("xref: "): self.xrefs.append( line[line.index(":")+2:] )
      if line.startswith("synonym: "): self.synonyms.append( line[line.index(":")+2:] )

  def get_id(self): return self.id
  def get_name(self): return self.name
  def get_namespace(self): return self.namespace
  def get_definition(self): return self.definition
  def get_is_a(self): return self.is_a
  def get_is_a_goids(self): return [e.split()[0] for e in self.is_a]
  def get_alt_ids(self): return self.alt_ids
  def get_xrefs(self): return self.xrefs
  def get_synonyms(self): return self.synonyms
  def get_is_obsolete(self): return self.obsolete
