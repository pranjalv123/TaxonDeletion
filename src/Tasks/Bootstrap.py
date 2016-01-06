"""
Xylem - Phylogenetic Pipelines with MPI

Bootstrap.py contains functions for bootstrapping 
    
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
import random

class BootstrapGenes(xylem.Task):
    def setup(self, nreps):
        self.nreps = nreps
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def desc(self):
        return str(self.maxlen)
    def run(self):
        mats = self.input_data["alignments"]
        output = dendropy.DnaCharacterMatrix(taxon_namespace = mats[0].taxon_namespace)
        output.fill_taxa()
        for i in range(len(mats)):
            x = random.randint(0, len(mats) - 1)
            output.extend_sequences(mats[x])
        self.result = {"alignments": [output]}
        return self.result
    
