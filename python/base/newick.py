import os, re, string
import needlemanwunsch

# =============================================================================
def get_motif_array(ssrarray):
  if len(ssrarray) == 0: return ssrarray
  if type(ssrarray[0]) is unicode: return ssrarray
  return [e.std_motif + '/' + e.gene_feature for e in ssrarray]

# =============================================================================
def get_consensus_from_aln(alns, removegaps=True):
  cons = []
  a1, a2 = alns[0], alns[1]
  for i in range(len(a1)):
    if a1[i] == '-' or a2[i] == '-':
      if not removegaps: cons.append("-")
    else: cons.append(a1[i])
  return cons

# =============================================================================
class Node():
  def __init__(self, name="", dist=0, parent=0):
    self.name = name
    self.distance_to_parent = dist
    if parent: self.set_parent(parent)
    else: self.parent = parent
    self.properties = {}
    self.children = []
    #print "init Node. name:", name, "dist", dist, "parent:", parent

  def summed_distance_to(self, other):
    if self == other: return 0
    self_nodes = [self]
    while 1:
      if self_nodes[-1].parent == 0: break
      self_nodes.append(self_nodes[-1].parent)
    other_nodes = [other]
    while 1:
      if other_nodes[-1].parent == 0: break
      other_nodes.append(other_nodes[-1].parent)
    distance = 0
    for node in self_nodes:
      if node in other_nodes: break
      distance += 2* int(node.distance_to_parent)
    return distance

  def set_parent(self, node):
    self.parent = node
    node.properties = {}
    if not self in node.children: node.children.append(self)
    if len(node.children) == 2: 
      node.name = '(' + string.join([c.name for c in node.children], ",") + ')'
      node.children = sorted(node.children, key=lambda e: e.name)

  def ssrs(self):
    if not self.properties.has_key('ssrs') or (len(self.properties['ssrs']) == 0 and len(self.children) == 2):
      score, pointers, aln = needlemanwunsch.align(get_motif_array(self.children[0].ssrs()), get_motif_array(self.children[1].ssrs()), -1, 2, -10)
      #print "D1", get_motif_array(self.children[0].ssrs())
      #print "D2", get_motif_array(self.children[1].ssrs())
      #print "DC", aln
      self.properties['ssrs'] = get_consensus_from_aln(aln)
      #print "+ 2 +", len(self.properties['ssrs'])
    return self.properties['ssrs']


# =============================================================================
class Tree():
  def __init__(self, file):
    self.name = os.path.split(file)[1]
    self.leaves = {}
    self.temp_ancestral_nodes = []
    self.ancestral_nodes = []
    self.build_from_file(file)

  def add_leave_node(self, node, leafname):
    self.leaves[leafname] = node

  def add_ancestral_node(self, node):
    self.temp_ancestral_nodes.append(node)

  def get_last_ancestral_node(self):
    return self.temp_ancestral_nodes[-1]  

  def remove_last_ancestral_node(self):
    self.ancestral_nodes.append(self.temp_ancestral_nodes.pop(-1))
    if len(self.temp_ancestral_nodes) > 0: self.ancestral_nodes[-1].set_parent(self.get_last_ancestral_node())

  def get_root_node(self):
    parent = self.leaves.values()[0].parent
    while parent.parent: parent = parent.parent
    return parent

  def get_ancestral_node_of(self, node1, node2):
    lineage1 = [node1]
    while lineage1[-1].parent: lineage1.append(lineage1[-1].parent)
    anode = node2
    while 1:
      if not anode in lineage1: anode = anode.parent
      else: break
    return anode

  def build_from_file(self, file):
    tree = open(file).readline().strip()
    while not tree.startswith(";"):
      if tree.startswith("("): 
        self.add_ancestral_node(Node())
        tree = tree[1:]
      elif tree.startswith("):"):
        self.get_last_ancestral_node().distance_to_parent = re.match("\):(\d+)", tree).group(1)
        tree = tree[2+len(str(self.get_last_ancestral_node().distance_to_parent)):]
        self.remove_last_ancestral_node()
      elif tree.startswith(","): 
        tree = tree[1:]
      elif re.match("([A-Za-z]+):(\d+)", tree):
        name, length = re.match("([A-Za-z]+):(\d+)", tree).groups()
        self.add_leave_node(Node(name, length, self.get_last_ancestral_node()), name)
        tree = tree[1+len(name)+len(str(length)):]
