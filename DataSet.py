import dendropy
import numpy as np
import tempfile
import subprocess
import StringIO


def read_multiphylip(path, taxon_namespace):
    lines = [i for i in open(path).readlines() if i.strip()]
    ind = 0
    matrices = []
    while ind < len(lines):
        # parse first line
        n = int(lines[0].split()[0])
        # read first n lines
        matrices.append(dendropy.DnaCharacterMatrix.get_from_string('\n'.join(lines[ind:ind+n+1]), 'phylip', taxon_namespace=taxon_namespace))
        ind += n + 1
    return matrices


class DataSet:
    def __init__(self, namespace, genetrees=None, seqs=None, speciestree=None):
        self.taxon_namespace = namespace
        self.genetrees = genetrees
        self.seqs = seqs
        self.speciestree = speciestree
    
    
        
    def delete_taxa_uniform(self, ndelete, genetrees=None, seqs=None, speciestree=None, maxlen=None):
        dsnew, deletion_list = self.delete_taxa(ndelete, genetrees, seqs, speciestree, maxlen, True)
        if self.speciestree and speciestree:
            newspeciestree = self.speciestree.clone()
            dsnew.speciestree=newspeciestree
            newspeciestree.prune_taxa(deletion_list)
        return dsnew

    def delete_taxa(self, ndelete, sigma=0, genetrees=None, seqs=None, maxlen=None, const=False):
        print sorted([int(i.label) for i in self.seqs[0]._taxon_sequence_map.keys()])
        print self.seqs[0].as_string('phylip')

        deletion_lists = []
        if const:
            dl = (np.random.choice(list(self.taxon_namespace), size=ndelete, replace=False))
            print dl
            for i in self.seqs:
                deletion_lists.append(dl)
        else:
            for i in self.seqs:
                taxa = [j for j in i if len(i[j])]
                nd = min(ndelete + np.random.randn() * sigma,  len([1 for j in i if len(i[j])]))
                deletion_lists.append(np.random.choice(taxa, size=nd, replace=False))
        print deletion_lists
        newgenetrees=None
        newseqs=None
        newspeciestree=None
        if self.genetrees and genetrees:
            newgenetrees = dendropy.TreeList(taxon_namespace = self.taxon_namespace)
            for (gt, deletion_list) in zip(self.genetrees, deletion_lists):
                gtnew = gt.clone()
                gtnew.prune_taxa(deletion_list)
                newgenetrees.append(gtnew)
        if self.seqs and seqs:
            newseqs = []
            for (seq, deletion_list) in zip(self.seqs, deletion_lists):
                seqnew = seq.clone()

                print sorted([int(i.label) for i in seqnew._taxon_sequence_map.keys()])
                print seqnew.as_string('phylip')

                
                seqnew.discard_sequences(deletion_list)
                if maxlen:
                    for i in seqnew:
                        del seqnew[i][maxlen:]
                newseqs.append(seqnew)
        if const:
            return DataSet(self.taxon_namespace, newgenetrees, newseqs, newspeciestree), dl
        else:
            return DataSet(self.taxon_namespace, newgenetrees, newseqs, newspeciestree)
    def assert_genetrees_sequences_same(self):
        assert(len(self.genetrees) == len(self.seqs))
        for (g, s) in zip(self.genetrees, self.seqs):
            print sorted(list(set(s)))
            print sorted(list(set([i.taxon for i in g.leaf_nodes()])))
#            assert(len(set(s) - set([i.taxon for i in g.leaf_nodes()])) == 0)
            
    def assert_consistency(self):
        if self.genetrees and self.seqs:
            self.assert_genetrees_sequences_same()

    def est_trees_fasttree(self):
        
        sio = StringIO.StringIO()
        
        for i in self.seqs:
            i.write_to_stream(sio, schema="phylip", suppress_missing_taxa=True)
        
        proc = subprocess.Popen(['fasttree', '-nt', '-gtr', '-nopr', '-gamma', '-n', str(len(self.seqs))], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
#        cmd = FastTreeCommandLine('fasttree', input=seq_f.name, output=genetree_f.name)
#        cmd()
        trees, err = proc.communicate(sio.getvalue())
        print err
        print "Trees:"
        print trees
        self.genetrees = dendropy.TreeList.get_from_string(trees, 'newick')
        
    
    def write(self, prefix, genetrees=None, seqs=None, speciestree=None):
        self.assert_consistency()
        if self.genetrees and genetrees:
            f = prefix + "_gtrees"
            self.genetrees.write(path=f, schema="newick")
        if self.seqs and seqs:
            f = prefix + "_seqs"
            stream = open(f, 'w')
            for i in self.seqs:
                i.write_to_stream(stream, schema="phylip")
        if self.speciestree and speciestree:
            f = prefix + "_stree"
            self.speciestree.write(path=f, schema="newick")
