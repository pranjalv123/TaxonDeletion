import dendropy
import Task

class WriteSpeciesTree(Task.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
        self.cache = False
        self.local = True
    def inputs(self):
        return [("speciestree", dendropy.Tree)]
    def outputs(self):
        return []
    def run(self):
        stree = self.input_data["speciestree"]
        stree.write(path=self.path, schema="newick")
        self.result = {}
        return self.result


class WriteGeneTrees(Task.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
        self.cache = False
        self.local = True
    def inputs(self):
        return [("genetrees", dendropy.Tree)]
    def outputs(self):
        return []
    def run(self):
        self.input_data["genetrees"].write(path=self.path, schema="newick")
        self.result = {}
        return self.result


class WritePhylip(Task.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
        self.cache = False
        self.local = True
    def inputs(self):
        return [("genetrees", dendropy.Tree)]
    def outputs(self):
        return []
    def run(self):
        stream = open(self.path, 'w')
        for al in self.input_data["alignments"]:
            al.write_to_stream(stream, schema="newick")
        self.result = {}
        return self.result
