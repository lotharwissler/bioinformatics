#!/usr/bin/python
import os, sys, re, string

class BibtexEntry:


  def __init__(self, lines):
    self.ATTRIBUTE_REGEX = re.compile("\s{2}(\S+)\s{1}=\s\{(.*)\}*$")
    self.BIBTEXSTART_REGEX = re.compile("@([A-Z]+)\{(\S+),$")
    self.key = ""
    self.bibtype = ""
    self.attributehash = {}
    while 1:
      if len(lines) == 0: break
      line = lines.pop(0)

      # end of entry
      if line.startswith("}"): break

      # bibtex entry start line and key definition
      if self.BIBTEXSTART_REGEX.match(line):
        self.bibtype = self.BIBTEXSTART_REGEX.match(line).group(1)
        self.key     = self.BIBTEXSTART_REGEX.match(line).group(2)
        continue

      # bibtex attribute start
      if self.ATTRIBUTE_REGEX.match(line):
        attr  = self.ATTRIBUTE_REGEX.match(line).group(1)
        value = self.ATTRIBUTE_REGEX.match(line).group(2)
        self.attributehash[attr] = value
      else: self.attributehash[attr] += " " + line.strip()

    for attr, value in self.attributehash.iteritems():
      if value.endswith("}"): self.attributehash[attr] = value[:-1]
      elif value.endswith("},"): self.attributehash[attr] = value[:-2]

  def get_key(self): return self.key
  def get_first_author(self): return self.attributehash['author'].split(" and ")[0]
  def get_attr(self, name):
    if self.attributehash.has_key(name): return self.attributehash[name]
    return ""

  def get_author_count(self, return_str=0): 
    count = self.attributehash['author'].count(" and ") +1
    if return_str: return "%s" % count
    else: return count

  def annotate(self):
    self.attributehash['annotate'] = "(%s co-authors)" % self.get_author_count()
    self.attributehash['note'] = "(%s co-authors)" % self.get_author_count()

  def to_s(self, escape_title=1, annotate=0):
    print "@" + self.bibtype + "{" + self.key + ","
    all_attrs = self.attributehash.keys()
    for i in range(len(all_attrs)):
      attr = all_attrs[i]
      if i == len(all_attrs)-1: comma = ""
      else: comma = ","
      if attr == "title" and escape_title: print "  " + attr + " = \"{" + self.attributehash[attr] + "}\"" + comma
      else: print "  " + attr + " = {" + self.attributehash[attr] + "}" + comma
    print "}"





def usage():
  print >> sys.stderr, "usage: " + sys.argv[0] + " db.bib  [n=max-coauthors]"
  sys.exit(1)


def plausi():
  if len(sys.argv) != 3: usage()
  inFile, inAuthors = sys.argv[1:3]
  return inFile, int(inAuthors)


def main():
  inFile, inAuthors = plausi()
  fo = open(inFile)
  while 1:
    line = fo.readline().rstrip()
    if line.startswith("%"): continue
    if line.startswith("@comment"): break
    if line.startswith("@"):
      lines = []
      lines.append(line)
      while 1:
        line = fo.readline().rstrip()
        lines.append(line)
        if line.startswith("}"): break
      b = BibtexEntry(lines)
      count = b.get_author_count()
      if count > inAuthors: b.annotate()
      b.to_s()

  fo.close()


main()
