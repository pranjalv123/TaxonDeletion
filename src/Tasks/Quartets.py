"""
Xylem - Phylogenetic Pipelines with MPI

Quartets.py contains routines for generating and manipulating weighted
quartet sets.
    
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
import dendropy
import copy
import tempfile

class WeightedQuartetSet:
    def __init__(self, tn, infile=None, mode='newick'):
        if isinstance(tn, WeightedQuartetSet):
            self.taxon_namespace = tn.taxon_namespace
            self.tn = tn.tn
            self.d = copy.deepcopy(tn.d)
            self.indices = tn.indices
            
            return 
        self.taxon_namespace = tn
        self.tn = tn
        self.d = np.zeros([len(tn), len(tn), len(tn), len(tn)])
        self.indices = {}
        for i, tax in enumerate(tn):
            self.indices[tax] = i
        if infile:
            for line in open(infile).readlines():
                line = line.translate(None, '()').replace(',', ' ').replace('|', ' ').replace(';', ' ').replace(':', ' ')
                a,b,c,d,w = line.split()
                self[(a,b,c,d)] = w
                
                
    def __getitem__(self, (a, b, c, d)):
        if isinstance(a, int):
            return self.d[a, b, c, d]
        ix = self.indices
        return self.d[ix[a], ix[b], ix[c], ix[d]]
        
    def __setitem__(self, (a, b, c, d), w):
        if isinstance(a, int):
            self.d[a, b, c, d] = w
            self.d[a, b, d, c] = w
            self.d[b, a, c, d] = w
            self.d[b, a, d, c] = w
            
            self.d[c, d, a, b] = w
            self.d[d, c, a, b] = w
            self.d[c, d, b, a] = w
            self.d[d, c, b, a] = w
            return
        
        ix = self.indices
        self.d[ix[a], ix[b], ix[c], ix[d]] = w
        self.d[ix[a], ix[b], ix[d], ix[c]] = w
        self.d[ix[b], ix[a], ix[c], ix[d]] = w
        self.d[ix[b], ix[a], ix[d], ix[c]] = w
        
        self.d[ix[c], ix[d], ix[a], ix[b]] = w
        self.d[ix[d], ix[c], ix[a], ix[b]] = w
        self.d[ix[c], ix[d], ix[b], ix[a]] = w
        self.d[ix[d], ix[c], ix[b], ix[a]] = w

        
    def write(self, f, style, translate=False):
        if style == 'wqmc':
            fstring = "%s,%s|%s,%s:%f\n"
        elif style == 'newick':
            fstring = "((%s,%s),(%s,%s); %f\n"
        if isinstance(f, str):
            f = open(f, 'w')

        def write_fn(w0, w1, w2, (a, b, c, d)):
            f.write(fstring % (a.label, b.label, c.label, d.label, w0))
            return w0
        def write_fn_trans(w0, w1, w2, (a, b, c, d)):
            ix = self.indices

            f.write(fstring % (str(ix[a]), str(ix[b]), str(ix[c]), str(ix[d]), w0))
            return w0

        if not translate:
            self.transform(write_fn)
        else:
            self.transform(write_fn_trans)
            return dict(enumerate(self.tn))


    
    def transform(self, fn):
        for a in range(len(self.tn)):
            at = self.tn[a]
            for b in range(a+1, len(self.tn)):
                bt = self.tn[b]
                for c in range(b+1, len(self.tn)):
                    ct = self.tn[c]
                    for d in range(c+1, len(self.tn)):
                        w0 = self[(a, b, c, d)]
                        w1 = self[(a, c, b, d)]
                        w2 = self[(a, d, c, b)]
                        dt = self.tn[d]
                        
                        self[(a, b, c, d)] = fn(w0, w1, w2, (at, bt, ct, dt))
                        self[(a, c, b, d)] = fn(w1, w0, w2, (at, ct, bt, dt))
                        self[(a, d, c, b)] = fn(w2, w0, w1, (at, dt, ct, bt))

import os


class TransformQuartets(xylem.Task):
    def setup(self, fn):
        self.fn = fn

    def inputs(self):
        return [("quartets", WeightedQuartetSet)]

    def outputs(self):
        return [("quartets", WeightedQuartetSet)]

    def run(self):
        self.input_data["quartets"].transform(self.fn)
        self.result = {"quartets":self.input_data["quartets"]}
        return self.result

class SVDQuartetFrequencies(xylem.Task):
    def setup(self):
        pass
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("quartets", WeightedQuartetSet)]
    
    def eliminate_header(self,lines):
        ind = (i[0] for i in enumerate(lines) if i[1][:3] == "SVD").next()  + 4
        return lines[ind:]
    def eliminate_footer(self,lines):
        ind = len(lines) - (i[0] for i in enumerate(reversed(lines)) if i[1][:4] == "Time").next() - 2
        return lines[:ind]

    def parse_line(self, line, taxon_namespace):
        tokens = line.split()
        if len(tokens) == 7:
            tokens = tokens[1:]

        taxa = [taxon_namespace[int(i) - 1] for i in [tokens[0], tokens[1], tokens[3], tokens[4]]]
        weight = tokens[5]
        return taxa, weight


    def run(self):
        dna = self.input_data['alignments'][0].clone(depth=2)

        taxon_namespace = dna.taxon_namespace

        del_list = []
        for t in taxon_namespace:
            if len(dna[t]) == 0:
                del_list.append(t)

        for t in del_list:
            taxon_namespace.remove_taxon(t)
        
        taxon_namespace.sort(key=lambda x: x.label)

        f = tempfile.NamedTemporaryFile()
        
        dna.write_to_path(f.name + '.nex', 'nexus')

        import subprocess

        pp = os.path.dirname(__file__) + '/parse_paup.sh'
    
        out, err = subprocess.Popen(["bash", pp, f.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

        lines = self.eliminate_header(out.split('\n'))
        lines = self.eliminate_footer(lines)
    
        quartets = WeightedQuartetSet(taxon_namespace)

        for l in lines:
            (a,b,c,d), w = self.parse_line(l, taxon_namespace)
            quartets[(a,b,c,d)] = w
        self.result = {"quartets":quartets}
        return self.result


