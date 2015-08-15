import dendropy
import numpy as np

def delete_taxa(seqsubset):
    nseqs = len([i for i in seqsubset if '-' not in i]) - 1
    n0 = int(seqsubset[0].split()[0])
    n1 = int(seqsubset[0].split()[1])
    modseq = []
    modseq.append(str(nseqs) + '  ' + str(n1).replace('\n', ''))
    modseq.extend([i.replace('\n', '') for i in seqsubset[1:] if '-' not in i])
    return modseq
    
    
def isboundary(s):
    splits = s.split()
    if len(splits) > 1:
        return splits[1].isdigit()
    return False

import sys
if __name__ == '__main__':
    if len(sys.argv) < 1:
        print "delete_taxa.py input "
        exit()
    seq_inp = sys.argv[1]
    
    if seq_inp:
        seqs = [i for i in open(seq_inp).readlines() if i]

        seq_gene_starts = [i[0] for i in enumerate(seqs) if isboundary(i[1])]
        seq_gene_starts.append(len(seqs))


        
#    outseq = open(sys.argv[2], 'w')
    outseq = sys.stdout
    i = 0
    
    for (s1, s2) in zip(seq_gene_starts[:-1], seq_gene_starts[1:]):
        s1 = seq_gene_starts[i]
        s2 = seq_gene_starts[i+1]
        modseq = delete_taxa(seqs[s1:s2])
        if modseq == False:
            outseq.close()
            exit()
        outseq.write('\n'.join(modseq))
        outseq.write('\n')
        i += 1
                
        
    if seq_inp:
        outseq.close()
       
 
