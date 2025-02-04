import dendropy
import xylem.Task
import os
from xylem.Tasks import *

class WriteSpeciesTree(xylem.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
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


class WriteQuartets(xylem.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
        self.local = True
    def desc(self):
        return self.path
    def inputs(self):
        return [("quartets", Quartets.WeightedQuartetSet)]
    def outputs(self):
        return []
    def run(self):
        stream = open(self.path, 'w')
        q = self.input_data["quartets"]
        q.write(stream, 'newick')
        self.result = {}
        return self.result



import fcntl

class WriteScore(xylem.Task):
    def setup(self, location, desc):
        self.path = location
        self.description = desc
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

class WriteAttrs(xylem.Task):
    def setup(self, location, desc, attrs, consts = None):

        self.path = location
        self.description = desc
        self.attrs = []
        self.types = []
        for i in attrs:
            if type(i) is tuple:
                self.attrs.append(i[0])
                self.types.append(i[1])
            else:
                self.attrs.append(i)
                self.types.append(float)
        self.cache = False
        self._is_result = True
        f = open(location, 'a')
        fcntl.lockf(f, fcntl.LOCK_EX)
        if os.path.getsize(location) == 0:
            f.write('description' + ',' + ','.join(self.attrs) + '\n')
        f.close()
        if consts:
            self.consts = consts
        else:
            self.consts = {}
        
    def desc(self):
        return self.path + ':' + self.description
    def inputs(self):
        return [i for i in zip(self.attrs, self.types) if i[0] not in self.consts]
    def outputs(self):
        return []
    def run(self):
        print "HERE Running WriteAttrs"
        print self.consts
        print self.input_data
        self.consts.update(self.input_data)
        print self.consts
        s = self.description + ',' + ','.join((str(self.consts[i]) for i in self.attrs))
        
        f = open(self.path, 'a')
        fcntl.lockf(f, fcntl.LOCK_EX)
        f.write(s)
        f.write('\n')
        f.close()
        self.result = {}
        return self.result
