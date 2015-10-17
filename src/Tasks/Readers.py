"""
Xylem - Phylogenetic Pipelines with MPI

Readers.py contains routines that read data from disk.
    
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

class ReadSpeciesTree(Task):
    def setup(self, location):
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
    def setup(self, location):
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
    def setup(self, location):
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
