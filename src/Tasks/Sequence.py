"""
Xylem - Phylogenetic Pipelines with MPI

Sequence.py contains routines for manipulating sequence data
    
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
import dendropy


def delete_empty_sequences(mat):
    dellist = []
    for tax in mat:
        if all([i == '-' for i in mat[tax].values()]):
            dellist.append(tax)
    mat.discard_sequences(dellist)
    return mat

class Concatenate(xylem.Task):
    def setup(self):
        pass
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def run(self):
        mats = self.input_data["alignments"]
        output = dendropy.DnaCharacterMatrix(taxon_namespace = mats[0].taxon_namespace)
        output.fill_taxa()
        for mat in mats:
            mat.pack('-')
            output.extend_sequences(mat)
        self.result = {"alignments": [delete_empty_sequences(output)]}
        return self.result

class PseudoGenes(xylem.Task):
    def setup(self, pgLength):
        self.pgLength = pgLength
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def run(self):
        mats = self.input_data["alignments"]
        concat = dendropy.DnaCharacterMatrix(taxon_namespace = mats[0].taxon_namespace)
        concat.fill_taxa()
        for mat in mats:
            mat.pack('-')
            concat.extend_sequences(mat)
            
        output = []

        for i in range(0, len(concat[0]), self.pgLength):
            subset = concat.new_character_subset('this_subset', range(i, i+self.pgLength))
            output.append(delete_empty_sequences(concat.export_character_subset('this_subset')))
            del(concat.character_subsets['this_subset'])
            
        self.result = {"alignments": output}
        return self.result
    
