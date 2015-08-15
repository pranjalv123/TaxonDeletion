import dendropy
import numpy as np

def delete_taxa(seqsubset, n, maxlen):
    splits = [i.split() for i in seqsubset]
    splits = [i for i in splits if i]
    alltaxa = [i[0] for i in splits if i[1][0] in 'ACTG']
    if len(alltaxa) == 0:
        return False
    l0 = int(seqsubset[0].split()[1])
    labels = set(np.random.choice(alltaxa, size=n, replace=False))
    modseq = []
    if l0 > maxlen and maxlen > 0:
        l = maxlen 
    else:
        l = l0
    modseq.append(seqsubset[0].replace(str(l0), str(l)).replace("\n", ""))
#    l = l + 7 #7 characters for sequence label
    for i in splits[1:]:
        if i[0] in labels:
            modseq.append(i[0] + " " + i[1].replace('A', '-').replace('C', '-').replace('T', '-').replace('G', '-')[:l])
        else:
            modseq.append(i[0] + " " + i[1][:l])
    return modseq
    
    
def isboundary(s):
    splits = s.split()
    if len(splits) > 1:
        return splits[1].isdigit()
    return False

import sys
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print "delete_taxa.py input maxlen ndelete...\nInput should be an alignment file in PHYLIP format."
        exit()
    seq_inp = sys.argv[1]
    ndelete = []
    maxlen = int(sys.argv[2])
    for i in sys.argv[3:]:
        ndelete.append(int(i))
    
    ndelete.sort()
    print seq_inp
    print ndelete
    if seq_inp:
        seqs = [i for i in open(seq_inp).readlines() if i]

        seq_gene_starts = [i[0] for i in enumerate(seqs) if isboundary(i[1])]
        seq_gene_starts.append(len(seqs))


    for n in ndelete:

        outseq = open(seq_inp + "l%02dm%02d" % (maxlen, n), 'w')
        i = 0

        for (s1, s2) in zip(seq_gene_starts[:-1], seq_gene_starts[1:]):
            s1 = seq_gene_starts[i]
            s2 = seq_gene_starts[i+1]
            modseq = delete_taxa(seqs[s1:s2], n, maxlen)
            if modseq == False:
                outseq.close()
                exit()
            outseq.write('\n'.join(modseq))
            outseq.write('\n')
            i += 1

        if seq_inp:
            outseq.close()
        
