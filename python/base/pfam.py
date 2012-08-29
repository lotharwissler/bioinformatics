def read_hmmout(ifile, evalue=10, matchreq=0.0):
  hash = {}
  fo = open(ifile)
  for line in fo:
    cols = line.strip().split()
    if len(cols) == 16: 
      i = []
      i.append(line.index("\t"))
      if line.count(" ") > 0: i.append(line.index(" "))
      line = line[min(i):]
    pd = PfamDomain(line)
    if float(pd.get_attr('E-value')) > evalue: continue
    if matchreq > 0 and ((float(pd.get_attr('alignment_end'))-float(pd.get_attr('alignment_start')))/float(pd.get_attr('hmm_length'))) < matchreq: continue
    #print pd.get_attr('seq_id'), pd.get_attr('hmm_name')
    #print pd.get_attr('seq_id')
    if not hash.has_key(pd.get_attr('seq_id')): hash[pd.get_attr('seq_id')] = []
    hash[pd.get_attr('seq_id')].append(pd)
  fo.close()
  return hash
  

class PfamDomain():
  def __init__(self, line):
    self.attributes = ['seq_id', 'alignment_start', 'alignment_end', 'envelope_start', 'envelope_end', 'hmm_acc', 'hmm_name', 'type', 'hmm_start', 'hmm_end', 'hmm_length', 'bit_score', 'E-value', 'significance', 'clan']
    line = line.strip()
    self.values =  line.split()

  def get_attr(self, name):
    if not name in self.attributes: return ""
    return self.values[ self.attributes.index(name) ]
  
  def covers(self, position):
    position = int(position)
    if int(self.get_attr('alignment_start')) <= position and int(self.get_attr('alignment_end')) >= position:
      return True
    return False
