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

import xylem.Task
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid
import dendropy

class ReadSpeciesTree(xylem.Task):
    def setup(self, location, *args, **kwargs):
        self.path = location
        self.cache = False
        self.local = True
    def desc(self):
        return self.path
    
    def inputs(self):
        return []
    def outputs(self):
        return [("speciestree", dendropy.Tree)]
    def run(self):
        s = open(self.path).read()
        if s[-1] != ';':
            s += ';'
        tree = dendropy.Tree.get_from_string(s, "newick")
        print "read species tree with", len(tree.taxon_namespace), "taxa"
        self.result = {"speciestree": tree}
        return self.result
    
class ReadGeneTrees(xylem.Task):
    def setup(self, location, limit=-1):
        self.path = location
        self.limit = limit
        self.cache = False
    def desc(self):
        return self.path
    def inputs(self):
        return []
    def outputs(self):
        return [("genetrees", dendropy.TreeList)]
    def run(self):

        if self.limit > 0:
            trees = dendropy.TreeList()
            for line in open(self.path).readlines()[:self.limit]:
                trees.append(dendropy.Tree.get_from_string(line, 'newick'))
        else:
            trees = dendropy.TreeList.get_from_path(self.path, "newick")

        print "read", len(trees), "gene trees with", len(trees.taxon_namespace), "taxa from", self.path
        self.result= {"genetrees": trees}
        return self.result


class ReadPhylip(xylem.Task):
    def setup(self, location, limit=-1):
        self.path = location
        self.limit = limit
        self.cache = False
    def desc(self):
        return self.path
    def inputs(self):
        return []
    def outputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def run(self):
        taxon_namespace = dendropy.TaxonNamespace()
        lines = [i for i in open(self.path).readlines() if i.strip()]
        ind = 0
        matrices = []
        if self.limit == -1:
            self.limit = len(lines)
        while ind < len(lines) and len(matrices) < self.limit:
            # parse first line
            n = int(lines[0].split()[0])
            # read first n lines
            matrices.append(dendropy.DnaCharacterMatrix.get_from_string('\n'.join(lines[ind:ind+n+1]), 'phylip', taxon_namespace=taxon_namespace))
            ind += n + 1
            # if len(matrices) % 10 == 0:
            #     print len(matrices), "sequences read"
        print "read", len(matrices), "sequences with", len(taxon_namespace), "taxa from", self.path
        self.result = {"alignments":matrices}
        return self.result

