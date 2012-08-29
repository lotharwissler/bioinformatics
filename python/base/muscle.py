import tempfile, os, fasta

# =============================================================================
def align(sequences, ids, outfile=False):
  h, infile = tempfile.mkstemp()
  os.close(h)
  fw = open(infile, 'w')
  for i in range(len(sequences)): fw.write(">" + ids[i] + "\n" + sequences[i] + "\n")
  fw.close()
  h, outfile = tempfile.mkstemp()
  os.close(h)
  os.system("muscle -in %s -out %s -quiet 2> /dev/null" %(infile, outfile))
  os.unlink(infile)
  aligned_sequences = []
  alnhash = fasta.get_sequence_hash(outfile)
  for gid in ids: aligned_sequences.append(alnhash[gid])
  os.unlink(outfile)
  return aligned_sequences
