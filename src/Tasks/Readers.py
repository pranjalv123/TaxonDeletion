import Task
import DataSet
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid

class ReadSpeciesTree(Task):
    def __init__(self, location, *args, **kwargs):
        super(ReadSpeciesTree, self).__init__(*args, **kwargs)
        self.path = location
        self.cache = False
        self.local = True
    def inputs(self):
        return []
    def outputs(self):
        return [("speciestree", dendropy.Tree)]
    def run(self):
        tree = dendropy.Tree.get_from_path(self.path, "newick")
        print "read species tree with", len(tree.taxon_namespace), "taxa"
        self.result = {"speciestree": tree}
        return self.result
    
class ReadGeneTrees(Task):
    def __init__(self, location, *args, **kwargs):
        super(ReadGeneTrees, self).__init__(*args, **kwargs)
        self.path = location
        self.cache = False
    def inputs(self):
        return []
    def outputs(self):
        return [("genetrees", dendropy.TreeList)]
    def run(self):
        trees = dendropy.TreeList.get_from_path(self.path, "newick")
        print "read", len(trees), "gene trees with", len(trees.taxon_namespace), "taxa"
        self.result= {"genetrees": trees}
        return self.result


class ReadPhylip(Task):
    def __init__(self, location, *args, **kwargs):
        super(ReadPhylip, self).__init__(*args, **kwargs)
        self.path = location
        self.cache = False
    def inputs(self):
        return []
    def outputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def run(self):
        taxon_namespace = dendropy.TaxonNamespace()
        lines = [i for i in open(self.path).readlines() if i.strip()]
        ind = 0
        matrices = []
        while ind < len(lines):
            # parse first line
            n = int(lines[0].split()[0])
            # read first n lines
            matrices.append(dendropy.DnaCharacterMatrix.get_from_string('\n'.join(lines[ind:ind+n+1]), 'phylip', taxon_namespace=taxon_namespace))
            ind += n + 1
            # if len(matrices) % 10 == 0:
            #     print len(matrices), "sequences read"
        print "read", len(matrices), "sequences with", len(taxon_namespace), "taxa"
        self.result = {"alignments":matrices}
        return self.result
