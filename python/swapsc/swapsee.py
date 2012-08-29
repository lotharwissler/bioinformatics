#!/usr/bin/python

import os, sys 				# low level handling, such as command line stuff
import string					# string methods available
import re							# regular expressions
import getopt					# comand line argument handling
import copy           # clone an object
from low import *			# custom functions, written by myself

DEBUG = 1

# =============================================================================	
def show_help( ):
  """ displays the program parameter list and usage information """
  stdout( "reads SWAPSC output and offers evaluation and plotting abilities.\n" )
  stdout( "usage: " + sys.argv[0] + " -c <path> -o <path>" )
  stdout( " " )
  stdout( " option    description" )
  stdout( " -h        help (this text here)" )
  stdout( " -f        swapsc output file" )
  stdout( " -g        group water and terrestrial branches" )
  stdout( " -r        remove overlaps" )
  stdout( " -m        modes to consider [default: \"PS,AdN,NS\"]" )

  stdout( " " )

  sys.exit(1)

# =============================================================================
def handle_arguments():
  """ verifies the presence of all necessary arguments and returns the processed set of arguments """
  if len ( sys.argv ) == 1:
    stderr( "no arguments provided." )
    show_help()	

  try: # check for the right arguments
    keys, values = getopt.getopt( sys.argv[1:], "hf:grm:" )
  except getopt.GetoptError:
    stderr( "invalid arguments provided." )
    show_help()

  args = { 'rmoverlaps': 1, 'modes': ['PS','AdN','NS'] }
  for key, value in keys:
    if key == '-f': args['swapscout'] = value
    if key == '-g': args['group'] = 1
    if key == '-r': args['rmoverlaps'] = 1
    if key == '-m': args['modes'] = value.split(',')
  
  if not args.has_key('swapscout'):
    stderr( "swapsc out file missing." )
    show_help()
  if not file_exists( args.get('swapscout') ):
    stderr( "swapsc out file does not exist." )
    show_help()

  return args

# =============================================================================
def resolve_overlap(con1, con2):
  priority = { 'PS':10, 'NS':10, 'AdN':8, 'HS':7, 'S':6, 'AdN + S':7 }
  # check if con2 is within con1
  if con1.mode == con2.mode:
    start = min([con1.start, con2.start])
    stop = max([con1.stop, con2.stop])
    con1.start = start
    con1.stop = stop
    return [con1]
  if con2.start > con1.start and con2.stop < con1.stop:
    if priority[ con2.mode ] < priority[ con1.mode ]: return [con1]
    elif priority[ con2.mode ] > priority[ con1.mode ]: 
      # break up in three
      con3 = copy.copy(con1)
      con4 = copy.copy(con1)
      con3.stop = con2.start - 1
      con4.stsart = con2.stop + 1
      return [con3,con2,con4]
    else: return [con1,con2]
  else: # con2 just overlaps with con1
    if priority[ con1.mode ] > priority[ con2.mode ]:
      con2.start = 1 + con1.stop
      if con2.start >= con2.stop: return [con1]
      else: return [con1, con2]
    elif priority[ con1.mode ] < priority[ con2.mode ]:
      con1.stop = con2.start - 1
      if con1.stop <= con1.start: return [con2]
      else: return [con1, con2]
    else:
      overlapfrom, overlapto = con2.start, con1.stop
      con1.stop = overlapfrom -1
      con2.start = overlapto +1
      if con1.stop <= con1.start: return [con2]
      elif con2.start >= con2.stop: return [con1]
      else: return [con1,con2]

# =============================================================================
# =============================================================================
class Swapscout:
  """
  this class can parse SWAPSC output and stores the information the following way:
  nseq: number of input sequences
  seqlength:          length (nt) of the input MSA
  seqhash:            seqnumber => nucleotidesequence
  nodehash:           for the ancestral sequences, their seqnumber => array of child seqnumbers
  parameterhash:      stores values for parameters estimated for the simulated data
  constrainthash:     stores all signals reported by SWAPSC
                      branch => region => { mode, p }
  summaryhash:        summarizing stats at the end of the SWAPSC out such as percentage of codons
                      for each mode of selection, and P(neutral sites)
  """
  def __init__(self, args):
    self.nseq = 0
    self.seqlength = 0
    self.seqhash = {}
    self.nodehash = {}
    self.parameterhash = {}
    self.constrainthash = {}
    self.summaryhash = {}
    self.args = args

  def parse(self, file):
    fo = open(file)
    lines = fo.readlines()
    fo.close()
    # nseq + seqlength
    line = lines.pop(0)
    self.nseq = int(re.match('Number of sequences =\s*(\d+)', line).group(1))
#    if DEBUG: print "number seq:" , self.nseq
    line = lines.pop(0)
    self.seqlength = int(re.match('Length of alignment \(nucleotides\) =\s*(\d+)', line).group(1))


#    if DEBUG: print "seq legnth:" , self.seqlength
    # input sequences
    line = lines.pop(0)
    line = lines.pop(0)
    count = 1
    while not re.match('$', line):
      name = line.rstrip()
      seq = lines.pop(0).rstrip()
      self.seqhash[count] = { 'name': name, 'seq': seq }
      count += 1
      line = lines.pop(0)
#    if DEBUG: print "seqhash:", self.seqhash

    # tree and branches
    while not re.match('Branches:', line):
      line = lines.pop(0)
    line = lines.pop(0)
    while not re.match('$', line):
      if re.match('\d+\s+:\s+\d+\.\.\.\d+', line):
        first = int( re.match('(\d+)', line).group(1) )
        second = int( re.search('\s+(\d+)\.\.\.', line).group(1) )
        third = int( re.search('\.\.\.(\d+)', line).group(1) )
        self.nodehash[first] = [second, third]
      line = lines.pop(0)
#    if DEBUG: print "nodehash:", self.nodehash

    # ancestral sequences
    while not re.match('Ancestral sequences inferred by MP:', line):
      line = lines.pop(0)
    line = lines.pop(0)
    line = lines.pop(0)
    while re.match('node', line):
      line = line.rstrip()
      elements = line.split()
      count = int( re.search('node(\d+):', elements[0]).group(1) )
      self.seqhash[count] = { 'seq': elements[1] } 
      line = lines.pop(0)
#    if DEBUG: print "seqhash:", self.seqhash

    # parameter estimates
    while not re.match('Parameter estimates using simulated data:', line):
      line = lines.pop(0)
    line = lines.pop(0)
    line = lines.pop(0)
    while not re.match('Numbers of the species follow the input', line):
      if re.match('$', line): 
        line = lines.pop(0)
        continue
      elements = line.split(';')
      for e in elements:
        key, value = re.search('(.*)\s+=\s+(\S+)', e).groups()
        self.parameterhash[key] = value
      line = lines.pop(0)
#    if DEBUG: print "parameters:", self.parameterhash

    # selective constraints
    while not re.match('=================', line):
      line = lines.pop(0)
    line = lines.pop(0)
    branch = ''
    while not re.match('\S+', line):
      if not re.search('\S+', line): 
        line = lines.pop(0)
        continue
      # branch definition
      if re.match("\s+\d+\.\.\d+$", line):
        branch = line.strip()
        self.constrainthash[branch] = {}
      else:
        col = line.rstrip().split()
        if string.join(col,'') != '-------':
          p = string.join(col[6:9])
          mode = col[9]
          if len(col) == 12:
            mode = string.join(col[9:])
          self.constrainthash[branch][col[0]] = { 'mode': mode, 'p':p }
      line = lines.pop(0)
#    if DEBUG: print "constrainthash:", self.constrainthash

    # summary
    while not re.match('Selective constraints', line):
      line = lines.pop(0)
    line = lines.pop(0)
    line = lines.pop(0)
    while re.search('\S', line):
      line = line.rstrip()
      col = line.split()
      self.summaryhash[col[0]] = {'% codons': col[1], 'mean Ka': col[2], 'mean Ks': col[3], 'mean W': col[4]}
      line = lines.pop(0)

    line = lines.pop(0)
    self.summaryhash['P(neutral sites)'] = re.search('(\S+)$', line).group(1)
#    if DEBUG: print "summary:", self.summaryhash

  def create_cluster(self):
      cluster = Cluster(self.args, self.nseq, self.seqlength)
      ids = self.seqhash.keys()
      ids.sort()
      for id in ids:
        hash = self.seqhash[id]
        node = Node(id)
        node.seq = hash['seq']
        if hash.has_key('name'): node.name = hash['name']
        if self.nodehash.has_key(id): node.children = self.nodehash[id]
        cluster.addNode(node)

      pnodes = self.nodehash.keys()
      pnodes.sort()
      i = 0
      while i < (len(pnodes)-1):
        pnode = int(pnodes[i])
        children = self.nodehash.get(pnode)
        for child in children:
          cild = int(child)
          cluster.branches['%s..%s' %(pnode,child)] = 1
          #sys.stderr.write("added branch: %s..%s\n" %(pnode,child))
        i += 1

      for branch, hash in self.constrainthash.iteritems():
        for region, details in hash.iteritems():
          start, stop = region.split('..')
          c = Constraint(branch, start, stop)
          c.mode = details['mode']
          c.p = details['p']
          cluster.addConstraint(c)

      return cluster

# =============================================================================
class Node:
  def __init__(self, id):
    self.id = int(id)
    self.name = ''
    self.seq = ''
    self.children = []

# =============================================================================
class Branch:
  def __init__(self, id, nodefrom, nodeto):
    self.id = id
    self.nodefrom = nodefrom
    self.nodeto = nodeto

# =============================================================================
class Constraint:
  def __init__(self, branch, start, stop):
    self.branch = branch
    self.start = int(start)
    self.stop = int(stop)
    self.mode = ''
    self.p = ''

  def __cmp__(self, other):
    return cmp(self.start, other.start)

  def is_unsignificant(self):
    if self.p == 'P > 0.05': return 1
    else: return 0

  def overlaps_with(self, con, debug=0):
    #if self.start < con.start and self.stop > con.start: return 1
    #else: return 0
    if debug:
      sys.stderr.write("overlap?\t%s\t%s\t%s\t%s\t%s\t%s" %(self.mode,self.start,self.stop,con.mode,con.start,con.stop))
    if (self.start < con.start and self.stop < con.start) or (self.start > con.stop and self.stop > con.stop): 
      if debug: sys.stderr.write("\tno\n")
      return 0
    else: 
      if debug: sys.stderr.write("\tyes\n")
      return 1

  def to_s(self):
    return string.join([self.branch,'%s' % self.start, '%s' % self.stop,self.mode,self.p],"\t")

# =============================================================================
class Gap:
  def __init__(self,branch,start,stop):
    self.branch = branch
    self.start = start
    self.stop = stop
  
  def __cmp__(self,other):
    return cmp(self.start,other.start)

  def overlaps_with(self,other):
    if (other.start < self.start and other.stop < self.start) or (other.start > self.stop and other.stop > self.stop): return 0
    else: return 1

# =============================================================================
# =============================================================================
# PANEL: width
#   TRACK: name
#     FEAT: range, label, color, glyph
#
class Cluster:
  """ store one complete SWAPSC output """
  def __init__(self, args, nseq, seqlength):
    self.args = args
    self.nseq = int(nseq)
    self.seqlength = int(seqlength)
    self.nodes = {}
    self.branches = {}
    self.constraints = {}
    self.gaps = {}

  def get_node_by_name(self,name):
    for node in self.nodes.values():
      if node.name == name: return node
    return None

  def addNode(self, node):
    if node.name == '' and node.children != []:
      node.name = '(%s, %s)' %( self.nodes.get(node.children[0]).name, self.nodes.get(node.children[1]).name)
    self.nodes[node.id] = node

  def addConstraint(self, c):
    if not self.constraints.has_key(c.branch): self.constraints[c.branch] = {}
    self.constraints[c.branch][c.start] = c

  def addGap(self,c):
    if not self.gaps.has_key(c.branch): self.gaps[c.branch] = {}
    self.gaps[c.branch][c.start] = c

  def reduce_to_significant_constraint(self):
    """ remove all with P > 0.05 or with unwanted modes """
    duplicate = {}
    for branch, hash in self.constraints.iteritems():
      for start, con in hash.iteritems():
        if not con.is_unsignificant() and con.mode in self.args.get('modes'):
          if not duplicate.has_key(branch): duplicate[branch] = {}
          duplicate[branch][start] = con
    self.constraints = duplicate

  def remove_overlapping_constraints(self):
    """ 
    - put together overlapping signals of identical mode
    - resolve overlaps according to priorities of different modes
    """
    NRconstraints = {}
    for branch, hash in self.constraints.iteritems():
      sortedstarts = hash.keys()
      sortedstarts.sort()
      sortedcons = []
      for key in sortedstarts: sortedcons.append(hash[key])
      pos = 0
      exit = 0
      while not exit:
        sortedcons.sort()
        if pos+1 == len(sortedcons):
          exit = 1
          continue
        firstcon = sortedcons[pos]
        nextcon = sortedcons[pos+1]
        if firstcon.overlaps_with(nextcon):
          #print "OVERLAP | 1: %s %s %s | 2: %s %s %s" %(firstcon.mode,firstcon.start,firstcon.stop,nextcon.mode,nextcon.start,nextcon.stop)
          newcons = resolve_overlap(firstcon,nextcon)
          #tmp = ""
          #for n in newcons:
          #  tmp += ' %s %s %s' %(n.mode, n.start, n.stop)
          #print "resolved into:", tmp
          sortedcons.pop(pos)
          sortedcons.pop(pos)
          newcons.sort()
          newcons.reverse()
          for c in newcons:
            sortedcons.insert(pos, c)
          pos = 0
        else:
          pos += 1
      self.constraints[branch] = {}
      for con in sortedcons:
        self.constraints[con.branch][con.start] = con

  def add_gaps(self):
    for branch in self.branches.keys():
      node1 = self.nodes.get( int( re.match('(\d+)\.\.', branch).group(1) ) )
      node2 = self.nodes.get( int( re.search('\.\.(\d+)$', branch).group(1) ) )
      # walk along the sequence, check if gap present in one of the sequences
      # if gap, create new or prolong already existing constraint
      gapstart, gapstop = None, None
      i = 0
      while i < self.seqlength:
        if node1.seq[i] == '-' or node2.seq[i] == '-':
          if not gapstart: gapstart = i
          gapstop = i
        else:
          if gapstart != None and gapstop != None:
            gap = Gap(branch, gapstart, gapstop)
            self.addGap(gap)
            gapstart, gapstop = None, None
        i += 1
      if gapstart != None and gapstop != None:
        gap = Gap(branch, gapstart, gapstop)
        self.addGap(gap)

  def get_ancestor_of(self,nodeid):
    for b in self.branches.keys():
      node1 = int( re.match("(\d+)\.\.\d+",b).group(1) )
      node2 = int( re.match("\d+\.\.(\d+)",b).group(1) )
      if nodeid == node2:
        return node1
    return None

  def group_branches(self):
    allbranches = self.branches.keys()
#    for branch in allbranches:
#      for g in self.gaps[branch].values():
#        sys.stderr.write("gap:\t%s\t%s\t%s" %(g.branch,g.start,g.stop) +"\n")
#    sys.stderr.write("\n")
 
    # water branches to group
    waterbranches = []
    pooc = self.get_node_by_name("P.oceanica")
    zoma = self.get_node_by_name("Z.marina")
    waterbranches.append('%s..%s' %( self.get_ancestor_of(zoma.id), zoma.id))
    waterbranches.append('%s..%s' %( self.get_ancestor_of(pooc.id), pooc.id))
    waterbranches.append('%s..%s' %( self.get_ancestor_of(self.get_ancestor_of(pooc.id)), self.get_ancestor_of(pooc.id)))
    #sys.stderr.write("water branches to group:\n")
    #for b in waterbranches: sys.stderr.write("  %s\n" % b)

    # terrestrial branches to group
    terrestrialbranches = []
    sorghum = self.get_node_by_name("S.bicolor")
    oryza = self.get_node_by_name("O.sativa")
    populus = self.get_node_by_name("P.trichocarpa")
    arath = self.get_node_by_name("A.thaliana")
    terrestrialbranches.append('%s..%s' %(self.get_ancestor_of(arath.id), arath.id))
    terrestrialbranches.append('%s..%s' %(self.get_ancestor_of(populus.id), populus.id))
    terrestrialbranches.append('%s..%s' %(self.get_ancestor_of(sorghum.id), sorghum.id))
    terrestrialbranches.append('%s..%s' %(self.get_ancestor_of(self.get_ancestor_of(sorghum.id)), self.get_ancestor_of(sorghum.id)))
    #sys.stderr.write("terrestrial branches to group:\n")
    #for b in terrestrialbranches: sys.stderr.write("  %s\n" % b)

    groups = [waterbranches,terrestrialbranches]
    newconstraints = {}
    newgaps = {}
    newbranches = {}
    i = 0
    while i < len(groups):
      g = string.join(groups[i], ' + ')
      newbranches[g] = 1
      newconstraints[g] = {}
      newgaps[g] = {}
      for branch in groups[i]:
        if not self.constraints.has_key(branch) and not self.gaps.has_key(branch): 
          #sys.stderr.write("branch %s is empty of cons and gaps\n" %(branch) )
          continue
        if self.constraints.has_key(branch) and not self.constraints[branch] == {}:
          cons = self.constraints[branch].values()
          for con in cons:
            con.branch = g
            add = 0
            if newconstraints[g].has_key(con.start):
              while newconstraints[g].has_key(con.start + add): add += 1
            newconstraints[g][con.start + add] = con
        if self.gaps.has_key(branch):
          gaps = self.gaps[branch].values()
          for gap in gaps:
            gap.branch = g
            add = 0
            if newgaps[g].has_key(gap.start):
              while newgaps[g].has_key(gap.start + add): add += 1
            newgaps[g][gap.start + add] = gap
      #print "GROUP", i, "|", g.replace('+',' + ')
      i += 1
    self.constraints = newconstraints
    self.branches = newbranches
    self.gaps = newgaps
#    for branch in self.branches.keys():
#      for g in self.gaps[branch].values():
#        sys.stderr.write("gap:\t%s\t%s\t%s" %(g.branch,g.start,g.stop) +"\n")
#    sys.stderr.write("\n")
                         
  def collapse_gaps(self):
    for branch in self.branches:
      gaps = self.gaps[branch].values()
      gaps.sort()
 #     for g in gaps:
 #       sys.stderr.write("gap:\t%s\t%s\t%s" %(g.branch,g.start,g.stop) +"\n")
      pos, exit = 0, 0
      while pos < (len(gaps) -1):
        gaps.sort()
        firstgap = gaps[pos]
        nextgap = gaps[pos+1]
        if firstgap.overlaps_with(nextgap):
          #sys.stderr.write(string.join([branch,str(firstgap.start),str(firstgap.stop),str(nextgap.start),str(nextgap.stop)],"\t"))
          # gap2 is within gap1 --> remove gap2
          if nextgap.start >= firstgap.start and nextgap.stop <= firstgap.stop: 
            gaps.pop(pos+1)
 #           sys.stderr.write(" | within --> remove gap2" +"\n")
          # gap2 overlaps with gap1 on the right side --> prolong gap1 and remove gap2
          else:
            firstgap.stop = nextgap.stop
            gaps.pop(pos+1)
 #           sys.stderr.write(" | overlap --> merging into new gap %s - %s" %(firstgap.start,firstgap.stop) +"\n")
          pos -= 1
          if pos < 0: pos = 0
        else:
          pos += 1
      self.gaps[branch] = {}
      for gap in gaps:
        self.gaps[branch][gap.start] = gap
  
  def remove_signals_in_gaps(self):
    newcons = {}
    gaps = []
    constraints = []
    for branch in self.branches: 
      gaps.extend( self.gaps[branch].values() )
      constraints.extend( self.constraints[branch].values() )
      newcons[branch] = {}
    
    for gap in gaps:
      pos = 0
      while pos < (len(constraints)):
        constraints.sort()
        con = constraints[pos]
        if gap.overlaps_with(con):
          #sys.stderr.write( string.join(["con",con.mode,str(con.start),str(con.stop),"gap",str(gap.start),str(gap.stop)],"\t") + "\n" )
   
          # con within gap --> remove con
          if gap.start <= con.start and gap.stop >= con.stop:
            constraints.pop(pos)
            #pos -= 1
            #if pos < 0: pos = 0
          else:
            #  gap within con --> split con
            if gap.start > con.start and gap.stop < con.stop:
              newcon = Constraint(con.branch, gap.stop+1, con.stop)
              newcon.mode = con.mode
              newcon.p = con.p
              con.stop = gap.start - 1
              constraints[pos] = con
              constraints.insert(pos,newcon)

            elif gap.start <= con.start:
              con.start = gap.stop +1
              constraints[pos] = con
            elif gap.start >= con.start:
              con.stop = gap.start -1
              constraints[pos] = con
        else:
          #sys.stderr.write("no overlap: GAP\t%s\t%s\tCON\t%s\t%s\n" % (gap.start, gap.stop, con.start,con.stop))
          pos += 1


    for con in constraints:
#      gap = Gap("",1023,1900)
#      if not gap.overlaps_with(con):
      if newcons[con.branch].has_key(con.start):
        oldcon = newcons[con.branch][con.start]
        if con.stop <= oldcon.stop: continue
      newcons[con.branch][con.start] = con

    self.constraints = newcons
  
  def reduce_to_unique_constraints(self):
    # build "consensus"
    seq = {}
    for branch in self.branches.keys():
      seq[branch] = self.seqlength*[""]
      for c in self.constraints[branch].values():
        for i in range(c.start-1,c.stop):
          if seq[branch][i] == "" or seq[branch][i] == c.mode: seq[branch][i] = c.mode
          else: seq[branch][i] = "X"
      #sys.stderr.write(branch + " | " + string.join(seq[branch], " ") + "\n\n" )

    # if any seq has X, everyone gets an X at pos i
    for i in range(self.seqlength):
      ambig = 0
      for seqs in seq.values():
        if seqs[i] == 'X': 
          ambig = 1
          break
      if ambig: 
        for branch in seq.keys(): 
          seq[branch][i] = 'X'

    # re-create constraints
    for branch in self.branches.keys():
      newconstraints = {}
      start, stop, mode = None, None, None
      for i in range(self.seqlength):
        s = seq[branch][i]
        if (s == "" or s == "X"):
          if start != None:
            c = Constraint(branch, start, stop)
            c.mode = mode
            newconstraints[c.start] = c
            start, stop, mode = None, None, None
          else: continue
        else:
          if start != None and s == mode:
            stop = i
          elif start != None and s != mode:
            c = Constraint(branch, start, stop)
            c.mode = mode
            newconstraints[c.start] = c
            start, stop, mode = i, i, s
          else:
            start, stop, mode = i, i, s
      if start != None and not newconstraints.has_key(start):
        c = Constraint(branch, start, stop)
        c.mode = mode
        newconstraints[c.start] = c
      self.constraints[branch] = newconstraints
    # TODO: ignore those parts that are ambiguous in one branch in the other branch as well

    # reduce to different signals
    branches = self.branches.keys()
    branch1 = branches[0]
    branch2 = branches[1]
    cons1 = self.constraints[branch1].values()
    cons2 = self.constraints[branch2].values()
    #sys.stderr.write("constraints in branch 1 (%s):\n" % branch1)
    #for con in cons1:
    #  sys.stderr.write("\t%s\t%s\t%s\n" %(con.mode,con.start,con.stop))
    #sys.stderr.write("constraints in branch 2 (%s):\n" % branch2)
    #for con in cons2:
    #  sys.stderr.write("\t%s\t%s\t%s\n" %(con.mode,con.start,con.stop))
    i = 0
    while i < (len(cons1)): 
      j = 0
      while j < (len(cons2)):
        if i < 0: i = 0
        if j < 0: j = 0
        if len(cons1) == 0: break
        #sys.stderr.write("i: %s j: %s ai: %s aj: %s\n" %(i,j, len(cons1), len(cons2)))
        con2 = cons2[j]
        con1 = cons1[i]
        if not con1.overlaps_with(con2):
          j += 1
          continue
        else:
          if con1.mode != con2.mode: 
            j += 1
            continue
          else:
            # con1 = con2
            if con1.start == con2.start and con1.stop == con2.stop:
              #sys.stderr.write("  1 == 2\n")
              cons1.pop(i)
              cons2.pop(j)

            # con2 within con1
            elif con2.start >= con1.start and con2.stop <= con1.stop:
              #sys.stderr.write(" 2 within 1\n")
              cons2.pop(j)
              cons1.pop(i)
              newstop = con2.start -1
              newstart = con2.stop +1
              con2.branch = con1.branch
              con2.start = newstart
              con2.stop = con1.stop
              con1.stop = newstop
              cons1.append(con1)
              cons1.append(con2)

            # con1 within con2
            elif con1.start >= con2.start and con1.stop <= con2.stop:
              #sys.stderr.write(" 1 within 2\n")
              cons2.pop(j)
              cons1.pop(i)
              newstart = con1.stop +1
              newstop = con1.start -1
              con1.branch = con2.branch
              con1.start = newstart
              con1.stop = con2.stop
              con2.stop = newstop
              cons2.append(con2)
              cons2.append(con1)

            # con1 left of con2
            elif con1.start <= con2.start:
              #sys.stderr.write("  con1 --> con2\n")
              #sys.stderr.write("   con1 init.: %s\n" % con1.to_s() )
              #sys.stderr.write("   con2 init.: %s\n" % con2.to_s() )
              newstart = con1.stop +1
              newstop = con2.start -1
              con1.stop = newstop
              con2.start = newstart
              cons1[i] = con1
              cons2[j] = con2
              #sys.stderr.write("   con1 modf.: %s\n" % con1.to_s() )
              #sys.stderr.write("   con2 modf.: %s\n" % con2.to_s() )

            # con1 right of con1
            elif con1.start >= con2.start:
              #sys.stderr.write("  con2 --> con1\n")
              #sys.stderr.write("   con1 init.: %s\n" % con1.to_s() )
              #sys.stderr.write("   con2 init.: %s\n" % con2.to_s() )

              newstart = con2.stop +1
              newstop = con1.start -1
              con2.stop = newstop
              con1.start = newstart
              cons1[i] = con1
              cons2[j] = con2
              #sys.stderr.write("   con1 modf.: %s\n" % con1.to_s() )
              #sys.stderr.write("   con2 modf.: %s\n" % con2.to_s() )


            i, j = -1, -1
      i += 1

    self.constraints = {}
    for con in cons1:
      if con.stop >= con.start:
        if not self.constraints.has_key(con.branch): self.constraints[con.branch] = {}
        self.constraints[con.branch][con.start] = con
    for con in cons2:
      if con.stop >= con.start:
        if not self.constraints.has_key(con.branch): self.constraints[con.branch] = {}
        self.constraints[con.branch][con.start] = con

    terrestrialcodons = 0
    aquaticcondons = 0
    aquatics = {'AdN':0, 'PS':0, 'NS':0}
    terrestrials = {'AdN':0, 'PS':0, 'NS':0}
    for branch in self.constraints.keys():
      for con in self.constraints[branch].values():
        add = int(round(1.0*(1 + con.stop - con.start)/3))
        if con.branch.count('+') == 2:
          aquaticcondons += add
          aquatics[con.mode] += add
        else:
          terrestrialcodons += add
          terrestrials[con.mode] += add
    # TODO: count gaps, then normalize codon count by the length of the region that was analyzable (non-gapped)
    gaps = []
    for branch in self.gaps.keys(): gaps.extend( self.gaps[branch].values() )
    tmpseq = [''] * self.seqlength
    for gap in gaps:
      for i in range(gap.start,gap.stop+1):
        tmpseq[i] = "g"
    ngaps = int(math.ceil(1.0*tmpseq.count("g") / 3))
    npcods = int(math.ceil(1.0*tmpseq.count("") / 3))

    cons = []
    for branch in self.constraints.keys(): cons.extend( self.constraints[branch].values() )
    tmpseq = [''] * self.seqlength
    for con in cons:
      for i in range(con.start,con.stop+1):
        tmpseq[i] = "c"
    ncods = int(math.ceil(1.0*tmpseq.count("c") / 3))
    #print self.seqlength / 3, (aquaticcondons + terrestrialcodons) /3, ngaps, npcods, ncods

    # output: name, total_length, analyzable_codons, constrained_codons, percent_constrained_codons, aquatic:AdN/PS/NS, terrestrial:AdN/PS/NS
    print "#" + string.join(['ID', 'aln_length', 'ungapped_length', 'constr_codons', 'constr_percent', 'aqua_codons', 'aqua_AdN', 'aqua_PS', 'aqua_NS', 'terr_codons', 'terr_AdN', 'terr_PS', 'terr_NS'], "\t")
    print "%s\t%s\t%s\t%s\t%01.2f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(self.args.get('swapscout'), self.seqlength/3, npcods, ncods, ncods*100.0/npcods, aquaticcondons, aquatics['AdN'], aquatics['PS'], aquatics['NS'], terrestrialcodons, terrestrials['AdN'], terrestrials['PS'], terrestrials['NS'])


  def debug(self):
    for branch, hash in self.constraints.iteritems():
      for start, con in hash.iteritems(): print con.to_s()

  def to_s(self):
    color = {}
    color['PS'] = "0.2,0.7,0.1"
    color['NS'] = "1,0.2,0.2"
    color['AdN'] = "0,0.5,1"
    color['HS'] = "0,0,0.7"
    color['S'] = "0.5,0.5,0.5"
    color['AdN + S'] = "1,0.5,0"
    color['gap'] = "0.6,0.6,0.6"

    flatfile = string.join(["PANEL", '%s' % self.seqlength],"\t") + "\n"
    # annotation
    annotationfilename = string.join(self.args.get('swapscout').split('.')[:3],'.') + '.ids.annotation'
    annotationline = 0
    if file_exists('./' + annotationfilename):
      annotationline = open(annotationfilename).readline().rstrip()
    if file_exists('../' + annotationfilename):
      annotationline = open('../' + annotationfilename).readline().rstrip()
    if annotationline:
      flatfile += string.join(["TRACK","gene description [TAIR8 annotation]","true"],"\t") + "\n"
      flatfile += string.join(["FEATURE",'%s..%s' %(1,self.seqlength -1),"1,0.8,0.2",annotationline.replace("\t"," ")],"\t") + "\n"
    
    # constraints and gaps
    for branch in self.branches.keys():
      if branch.find('+') == -1:
        node1 = int( re.match('(\d+)\.\.', branch).group(1) )
        node2 = int( re.search('\.\.(\d+)$', branch).group(1) )
        branchname = self.nodes[node1].name + ' : ' + self.nodes[node2].name
      else:
        if branch.count('+') >= 3: branchname = "terrestrial"
        else: branchname = "aquatic"
        #branchname = branch
      flatfile += string.join(["TRACK",branchname,"false"],"\t") + "\n"

      #if not self.constraints.has_key(branch): continue
      if self.constraints.has_key(branch):
        hash = self.constraints[branch]
        for start, con in hash.iteritems():
          flatfile += string.join(["FEATURE",'%s..%s' %(con.start,con.stop),color[con.mode]],"\t") + "\n"

      if self.gaps.has_key(branch):
        hash = self.gaps[branch]
        for start, gap in hash.iteritems():
          flatfile += string.join(["FEATURE",'%s..%s' %(gap.start,gap.stop),color['gap']],"\t") + "\n"

    # legend
    legend = {
      'PS': 'Positive selection',
      'HS': 'Hot spots',
      'S' : 'Saturation of synonymous sites',
      'AdN': 'Acceleration of non-synonymous substitutions',
      'NS': 'Negative selection',
      'AdN + S': 'Acceleration of non-syn. substitutions + Saturation of syn. sites'
    }
    flatfile += string.join(["TRACK","Legend","true"],"\t") + "\n"
    legendstartpos = int(0.01 * self.seqlength)
    legendstoppos = int(0.99 * self.seqlength)
    legendrange = '%s..%s' % (legendstartpos,legendstoppos)
    flatfile += string.join(["FEATURE", legendrange, color['gap'], "gap in at least one of the sequences"],"\t") + "\n"
    for m in self.args.get('modes'):
      flatfile += string.join(["FEATURE", legendrange, color[m], '%s: %s' %(m,legend[m])],"\t") + "\n"
    
    return flatfile


  def plot(self):
    count = 1
    outfile = self.args.get('swapscout') + '.%s.png' % count
    while file_exists(outfile):
      count += 1
      outfile = self.args.get('swapscout') + '.%s.png' % count
    output = self.to_s()
    tmpfile = ".bio-graphics-plot.txt"
    write_to_file(tmpfile,output)
    os.system("bio-graphics-plot.rb %s %s" %(tmpfile,outfile) )

          
# =============================================================================
# =============================================================================
def main( args ):
  so = Swapscout(args)
  file = args.get('swapscout')
  so.parse( file )
  cluster = so.create_cluster()
  #print 10*"=", "all signals", 10*"="
  #cluster.to_s()
  cluster.reduce_to_significant_constraint()
  #print 10*"=", "only significant", 10*"="
  #print cluster.to_s()
  cluster.remove_overlapping_constraints()
  #print 10*"=", "no overlaps", 10*"="
  #cluster.to_s()

  cluster.add_gaps()
  #print 10*"=", "with gaps", 10*"="
  #cluster.to_s()

  cluster.plot()
  #print 10*"=", "SWAPSEE", 10*"="
  if args.has_key('group'):
    cluster.group_branches()
    cluster.collapse_gaps()
    cluster.remove_signals_in_gaps()
    cluster.plot()
    cluster.reduce_to_unique_constraints()
    cluster.plot()


# =============================================================================
# === MAIN ====================================================================
# =============================================================================

args = handle_arguments(  )
main( args )
