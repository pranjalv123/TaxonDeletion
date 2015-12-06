import dendropy
import xylem.Task

class WriteSpeciesTree(xylem.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
        self.cache = False
        self.local = True
    def desc(self):
        return self.path
    def inputs(self):
        return [("speciestree", dendropy.Tree)]
    def outputs(self):
        return []
    def run(self):
        stree = self.input_data["speciestree"]
        stree.write(path=self.path, schema="newick")
        self.result = {}
        return self.result


class WriteGeneTrees(xylem.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
        self.cache = False
        self.local = True
    def desc(self):
        return self.path
    def inputs(self):
        return [("genetrees", dendropy.TreeList)]
    def outputs(self):
        return []
    def run(self):
        self.input_data["genetrees"].write(path=self.path, schema="newick")
        self.result = {}
        return self.result


class WritePhylip(xylem.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
        self.cache = False
        self.local = True
    def desc(self):
        return self.path
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return []
    def run(self):
        stream = open(self.path, 'w')
        for al in self.input_data["alignments"]:
            al.write_to_stream(stream, schema="phylip", suppress_missing_taxa=True)
        self.result = {}
        return self.result


import fcntl

class WriteScore(xylem.Task):
    def setup(self, location, desc):
        self.path = location
        self.description = desc
        self.cache = False
        self._is_result = True
        f = open(location, 'w')
        fcntl.lockf(f, fcntl.LOCK_UN)
    def desc(self):
        return self.path + ':' + self.description
    def inputs(self):
        return [("score", float)]
    def outputss(self):
        return []
    def run(self):
        f = open(self.path, 'a')
        fcntl.lockf(f, fcntl.LOCK_EX)
        f.write(self.description + ',' + self.input_data['score'])
        f.close()
        self.result = {}
        return self.result
