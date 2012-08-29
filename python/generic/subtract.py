#!/usr/bin/python

import sets
import sys, os

def get_lines( file ):
  lines = []
  fo = open(file)
  for line in fo:
    line = line.rstrip()
    lines.append(line)

  return sets.Set(lines)

ref = get_lines(sys.argv[1])
for filename in sys.argv[2:]:
  l = get_lines(filename)
  for e in l:
    if e in ref: ref.remove(e)

for e in ref: print e
