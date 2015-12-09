"""
Xylem - Phylogenetic Pipelines with MPI

Delete.py contains tasks that delete taxa from datasets.
    
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
import sys

class DeleteTaxaUniform(xylem.Task):
    def setup(self, ndelete, gtrees=True, alignments=True, stree=True, *args, **kwargs):
        self.ndelete = ndelete
        self.gtrees=gtrees
        self.alignments=alignments
        self.stree = stree
    def inputs(self):
        inp = []
        if self.gtrees:
            inp +=("genetrees", dendropy.TreeList)
        if self.seqs:
            inp +=  ("alignments", (dendropy.DnaCharacterMatrix,))
        if self.stree:
            inp += ("speciestree", dendropy.Tree)
        return inp
    def outputs(self):
        return self.inputs()
    def desc(self):
        return str(self.ndelete)
    def run(self):
        tn = None
        if self.seqs:
            dna = [i.clone() for i in self.input_data["alignments"]]
            tn = dna[0].taxon_namespace
        if self.gtrees:
            gt = self.input_data["genetrees"].clone()
            if tn:
                gt.migrate_taxon_namespace(tn)
            else:
                tn = gt.taxon_namespace
        if self.stree:
            st = self.input_data["speciestree"].clone()
            if tn:
                st.migrate_taxon_namespace(tn)
            else:
                tn = st.taxon_namespace

        deletion_list = np.random.choice(list(tn), size=self.ndelete, replace=False)
        self.result = []
        if self.gtrees:
            for g in gt:
                g.prune_taxa(deletion_list)
            self.result['genetrees'] = gt
        if self.seqs:
            for seq in dna:
                seq.discard_sequences(deletion_list)
            self.result['alignments'] = dna
        if self.stree:
            st.prune_taxa(deletion_list)
            self.result['speciestree'] = st
        return self.result
    
class DeleteTaxaRandom(xylem.Task):
    def setup(self, ndelete, sigma=0, *args, **kwargs):
        self.ndelete = ndelete
        self.sigma = sigma
    def desc(self):
        return str(self.ndelete) + ' +/- ' + str(self.sigma)
    def inputs(self):
        return [("genetrees", dendropy.TreeList), ("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("genetrees", dendropy.TreeList), ("alignments", (dendropy.DnaCharacterMatrix,))]
    def run(self):
        debug = False
        if '--debug' in sys.argv:
            debug = True
        dna = [i.clone() for i in self.input_data["alignments"]]
        gt = self.input_data["genetrees"].clone()
        gt.migrate_taxon_namespace(dna[0].taxon_namespace)
        for seq, g in zip(dna, gt):
            if debug:
                print g
            taxon_list = [i.taxon for i in g.leaf_nodes()]
                            
            nd = min(self.ndelete + np.random.randn() * self.sigma,  len(taxon_list) - 4)
            if debug:
                print nd
                print taxon_list
            deletion_list = np.random.choice(taxon_list, size=nd, replace=False)
            if debug:
                print deletion_list
            g.prune_taxa(deletion_list)
            if debug:
                print g
                print
            seq.discard_sequences(deletion_list)
        self.result = {"alignments":dna, "genetrees":gt}
        return self.result

class LimitSeqLength(xylem.Task):
    def setup(self, maxlen, *args, **kwargs):
        self.maxlen = maxlen
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def desc(self):
        return str(self.maxlen)
    def run(self):
        dna = [i.clone() for i in self.input_data["alignments"]]
        for seq in dna:
            for tax in seq:
                del seq[tax][self.maxlen:]
        self.result = {"alignments":dna}
        return self.result
    
