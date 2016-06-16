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
import xylem.Task
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid
import dendropy

class CompareTrees(xylem.Task):
    def setup(self, outputfile=None, tag="", *args, **kwargs):
        self.local=True
        self.outputfile = outputfile
        self.tag = tag
        self._is_result=True
    def desc(self):
        return ""
    def inputs(self):
        return[("truespeciestree", dendropy.Tree), ("estimatedspeciestree", dendropy.Tree)]
    def outputs(self):
        return[("rfdistance", float)]
    def run(self):
        truetree = self.input_data["truespeciestree"]
        estimatedtree = self.input_data["estimatedspeciestree"]
        tn = estimatedtree.taxon_namespace
        truetree.migrate_taxon_namespace(tn)
        truetree.update_bipartitions()
        estimatedtree.update_bipartitions()
        diff = dendropy.calculate.treecompare.false_positives_and_negatives(truetree, estimatedtree)[0]
        print len(truetree.encode_bipartitions())
        print len(estimatedtree.encode_bipartitions())
        print len(truetree.leaf_nodes())
        print len(estimatedtree.leaf_nodes())

        print dendropy.calculate.treecompare.false_positives_and_negatives(truetree, estimatedtree)
        print dendropy.calculate.treecompare.symmetric_difference(truetree, estimatedtree)
        print diff, len(truetree.internal_edges())
        print diff/float(len(truetree.internal_edges()))
        self.result ={"rfdistance":diff/float(len(truetree.internal_edges()))}
        if self.outputfile:
            open(self.outputfile, 'a').write(self.tag + ',' + str(diff) + '\n')
            
        return self.result


class CalculateAD(xylem.Task):
    def inputs(self):
        return [("truespeciestree", dendropy.Tree), ("genetrees", dendropy.TreeList)]
    def outputs(self):
        return [("AD", float)]
    def run(self):
        st = self.input_data["truespeciestree"]
        gt = self.input_data["genetrees"]
        
        distances = []

        for t in gt:
            st1 = st.clone()
            t1 = t.clone()
            t1.migrate_taxon_namespace(st1.taxon_namespace)
            print [i.taxon for i in t1.leaf_node_iter()]
            print [i.taxon for i in st1.leaf_node_iter()]
            st1.retain_taxa([i.taxon for i in t1.leaf_node_iter()])
            t1.retain_taxa([i.taxon for i in st1.leaf_node_iter()])

            distances.append(float(dendropy.calculate.treecompare.symmetric_difference(t1, st1))/(2.0 * len(t1.leaf_nodes())))

        self.result = {'AD':np.mean(distances)}
        return self.result
