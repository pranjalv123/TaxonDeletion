"""
Xylem - Phylogenetic Pipelines with MPI

Methods.py contains routines for interacting with phylogenetic methods.
    
Copyright (C) 2015 Pranjal Vachaspati
pr@nj.al

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import Task
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid
import ASTRID
import dendropy
import tempfile

class RunFastTree(Task.Task):
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("genetrees", dendropy.TreeList)]
    
    def desc(self):
        return ""
    
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


    
class RunASTRID(Task.Task):
    def setup(self, distmethod=None, *args, **kwargs):
        self.distmethod = distmethod
        if self.distmethod==None:
            self.distmethod="auto"
    def desc(self):
        return self.distmethod
    
    def inputs(self):
        return [("genetrees", dendropy.TreeList)]
    def outputs(self):
        return [("estimatedspeciestree", dendropy.Tree)]
    def run(self):
        a = ASTRID.ASTRID(self.input_data["genetrees"])
        a.run(self.distmethod)
        self.result = {"estimatedspeciestree": a.tree}
        print a.tree
        print type(a.tree)
        return self.result

class RunASTRAL(Task.Task):
    def setup(self, *args, **kwargs):
        pass
    def desc(self):
        return ""
    def inputs(self):
        return [("genetrees", dendropy.TreeList)]
    def outputs(self):
        return [("estimatedspeciestree", dendropy.Tree)]
    def run(self):

        print "RUNNING ASTRAL"
        
        f = tempfile.NamedTemporaryFile()

        gt = self.input_data["genetrees"]

        gt = dendropy.TreeList([i for i in gt if len(i.leaf_nodes()) > 3])
        
        gt.write(path=f.name, schema='newick')
        
        proc = subprocess.Popen(['ASTRAL', '-i', f.name], stdout=subprocess.PIPE)

        
        streestr, err = proc.communicate()
        print err
        print streestr
        stree = dendropy.Tree.get_from_string(streestr, 'newick')
        self.result = {"estimatedspeciestree": stree}
        return self.result
