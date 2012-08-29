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

l1 = get_lines(sys.argv[1])
l2 = get_lines(sys.argv[2])
for e in l1.difference(l2):
  print e
