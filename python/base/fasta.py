import re
import gzip 

# =============================================================================
def get_sequence_hash(fastafile):
  seqhash = {}
  key = ""
  if fastafile.endswith('.gz'): fo = gzip.open(fastafile)
  else: fo = open(fastafile)
  for line in fo:
    if line.startswith(">"):
      gid = re.match(">(\S+)", line).group(1)
      key = gid
      seqhash[key] = ""
    else:
      if key != "": seqhash[key] += line.strip()
  fo.close()
  return seqhash
  
# =============================================================================
def get_length_hash(fastafile):
  lenhash = {}
  key = ""
  if fastafile.endswith('.gz'): fo = gzip.open(fastafile)
  else: fo = open(fastafile)
  for line in fo:
    if line.startswith(">"):
      gid = re.match(">(\S+)", line).group(1)
      key = gid
      lenhash[key] = 0
    else:
      if key != "": lenhash[key] += len(line.strip())
  fo.close()
  return lenhash
