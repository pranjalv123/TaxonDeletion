import dendropy
import numpy as np

def delete_taxa(seqsubset, n, maxlen, filterfun):
    splits = [i.split() for i in seqsubset]
    splits = [i for i in splits if i]
    alltaxa = [i[0] for i in splits if i[1][0] in 'ACTG']
    if len(alltaxa) == 0:
        return False
    l0 = int(seqsubset[0].split()[1])

    labels = set(filterfun(alltaxa, n))
    modseq = []
    if maxlen >= 0 and  l0 > maxlen and maxlen > 0:
        l = maxlen 
    else:
        l = l0
    modseq.append(seqsubset[0].replace(str(l0), str(l)).replace("\n", ""))
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

funMap = {}

def randomfilter(alltaxa, n):
    labels = set(np.random.choice(alltaxa, size=n, replace=False))
    return labels

funMap['random'] = randomfilter


def constantfilter(alltaxa, n):
    if n not in constantfilter.labels:
        constantfilter.labels[n] = set(np.random.choice(alltaxa, size=n, replace=False))
    return constantfilter.labels[n]
constantfilter.labels = {}

funMap['constant'] = constantfilter

import sys
import argparse
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Delete taxa from PHYLIP sequence files")

    parser.add_argument('-i', '--input', nargs='*', action='append')
    parser.add_argument('-l', '--maxlen', dest='maxlen', nargs='*', action='append', type=int)
    parser.add_argument('-n', '--ndelete', dest='ndelete',nargs='*',  action='append', type=int)
    parser.add_argument('-t', '--type', dest='types',nargs='*', action='append')
    
    # if len(sys.argv) < 4:
    #     print "delete_taxa.py input maxlen ndelete...\nInput should be an alignment file in PHYLIP format."
    #     exit()
    # seq_inp = sys.argv[1]
    # ndelete = []
    # maxlen = int(sys.argv[2])
    # for i in sys.argv[3:]:
    #     ndelete.append(int(i))

    args = vars(parser.parse_args())

    ndelete = args['ndelete'][0]
    maxlens = args['maxlen'][0]
    if not maxlens:
        maxlens = [-1]

    inputs = args['input'][0]

    types = args['types'][0]
    
    ndelete.sort()
    print args
    print ndelete

    for seq_inp in inputs:
        if seq_inp:
            seqs = [i for i in open(seq_inp).readlines() if i]
            
            seq_gene_starts = [i[0] for i in enumerate(seqs) if isboundary(i[1])]
            seq_gene_starts.append(len(seqs))

        for maxlen in maxlens:
            for type_str in types:
                type = funMap[type_str]
                done = False
                for n in ndelete:
                    if done:
                        break
                    if maxlen >= 0:
                        outseq = open(seq_inp + "l%02dm%02d" % (maxlen, n), 'w')
                    else:
                        outseq = open(seq_inp + "m%02d" % (n, ), 'w')

                    i = 0

                    for (s1, s2) in zip(seq_gene_starts[:-1], seq_gene_starts[1:]):
                        s1 = seq_gene_starts[i]
                        s2 = seq_gene_starts[i+1]
                        modseq = delete_taxa(seqs[s1:s2], n, maxlen, type)
                        if modseq == False:
                            outseq.close()
                            done = True
                            break
                        outseq.write('\n'.join(modseq))
                        outseq.write('\n')
                        i += 1

        if seq_inp:
            outseq.close()

                    
