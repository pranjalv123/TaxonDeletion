import Task
import DataSet
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid
import ASTRID

class RunFastTree(Task):
    def __init__(self, *args, **kwargs):
        super(RunFastTree, self).__init__(*args, **kwargs)
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("genetrees", dendropy.TreeList)]
    def run(self):
        self.seqs = self.input_data["alignments"]
        sio = StringIO.StringIO()
        for seq in self.seqs:
            seq.write_to_stream(sio, schema="phylip", suppress_missing_taxa=True)
        proc = subprocess.Popen(['fasttree', '-nt', '-gtr', '-nopr', '-gamma', '-n', str(len(self.seqs))], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        trees, err = proc.communicate(sio.getvalue())
        print err
        genetrees = dendropy.TreeList.get_from_string(trees, 'newick')
        self.result = {"genetrees":genetrees}
        return self.result

class RunASTRID(Task):
    def __init__(self, distmethod=None, *args, **kwargs):
        super(RunASTRID, self).__init__(*args, **kwargs)
        self.distmethod = distmethod
        if self.distmethod==None:
            self.distmethod="auto"
    def inputs(self):
        return [("genetrees", dendropy.TreeList)]
    def outputs(self):
        return [("estimatedspeciestree", dendropy.Tree)]
    def run(self):
        a = ASTRID.ASTRID(self.input_data["genetrees"])
        a.run(self.distmethod)
        self.result = {"estimatedspeciestree": a.tree}
        return self.result
