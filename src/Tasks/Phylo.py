"""
Xylem - Phylogenetic Pipelines with MPI

Phylo.py contains tasks that handle common and simple phylogenetic
tasks.

    
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
import DataSet
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid

class CompareTrees(Task):
    def __init__(self, *args, **kwargs):
        super(CompareTrees, self).__init__(*args, **kwargs)
        self.local=True
    def inputs(self):
        return[("truespeciestree", dendropy.Tree), ("estimatedspeciestree", dendropy.Tree)]
    def outputs(self):
        return[("rfdistance", float)]
    def run(self):
        truetree = self.input_data["truespeciestree"]
        estimatedtree = self.input_data["estimatedspeciestree"]
        tn = estimatedtree.taxon_namespace
        truetree.migrate_taxon_namespace(tn)
        diff = dendropy.calculate.treecompare.false_positives_and_negatives(truetree, estimatedtree)[0]
        self.result ={"rfdistance":diff}
        print "RF distance is", diff
        return self.result
