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
import time

class WeightedQuartetSet:
    def __init__(self, tn, infile=None, mode='newick'):
        if isinstance(tn, WeightedQuartetSet):
            self.taxon_namespace = tn.taxon_namespace
            self.tn = tn.tn
            self.d = copy.deepcopy(tn.d)
            self.rejected = copy.deepcopy(tn.rejected)
            self.indices = tn.indices
            return 
        self.taxon_namespace = tn
        self.tn = tn
        self.d = np.zeros([len(tn), len(tn), len(tn), len(tn)])
        self.rejected = np.zeros([len(tn), len(tn), len(tn), len(tn)], dtype=np.bool_)
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
        
    def __setitem__(self, (a, b, c, d), w, rejected=False):
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

        if rejected:
            if isinstance(a, int):
                self.rejected[a, b, c, d] = True
                self.rejected[a, b, d, c] = True
                self.rejected[b, a, c, d] = True
                self.rejected[b, a, d, c] = True

                self.rejected[c, d, a, b] = True
                self.rejected[d, c, a, b] = True
                self.rejected[c, d, b, a] = True
                self.rejected[d, c, b, a] = True
                return

            ix = self.indices
            self.rejected[ix[a], ix[b], ix[c], ix[d]] = True
            self.rejected[ix[a], ix[b], ix[d], ix[c]] = True
            self.rejected[ix[b], ix[a], ix[c], ix[d]] = True
            self.rejected[ix[b], ix[a], ix[d], ix[c]] = True

            self.rejected[ix[c], ix[d], ix[a], ix[b]] = True
            self.rejected[ix[d], ix[c], ix[a], ix[b]] = True
            self.rejected[ix[c], ix[d], ix[b], ix[a]] = True
            self.rejected[ix[d], ix[c], ix[b], ix[a]] = True

            
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


    
    def transform(self, fn, reject=False):
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


                        rejected = self.rejected[a,b,c,d]
                        dt = self.tn[d]
                        
                        if reject:
                            self[(a, b, c, d)] = fn(w0, w1, w2, (at, bt, ct, dt), rejected)
                            self[(a, c, b, d)] = fn(w1, w0, w2, (at, ct, bt, dt), rejected)
                            self[(a, d, c, b)] = fn(w2, w0, w1, (at, dt, ct, bt), rejected)
                        else:
                            self[(a, b, c, d)] = fn(w0, w1, w2, (at, bt, ct, dt))
                            self[(a, c, b, d)] = fn(w1, w0, w2, (at, ct, bt, dt))
                            self[(a, d, c, b)] = fn(w2, w0, w1, (at, dt, ct, bt))

import os


class TransformQuartets(xylem.Task):
    def setup(self, fn, reject=False):
        self.fn = fn
        self.reject = reject

    def inputs(self):
        return [("quartets", WeightedQuartetSet)]

    def outputs(self):
        return [("quartets", WeightedQuartetSet)]

    def run(self):
        self.input_data["quartets"].transform(self.fn, self.reject)
        self.result = {"quartets":self.input_data["quartets"]}
        return self.result


class SVDQuartetSet(xylem.Task):
    def setup(self):
        pass
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]        
    def outputs(self):
        return [("quartets", WeightedQuartetSet)]

    def parse_line(self, line, taxon_namespace):
        halves = line.strip().split('|')
        tokens = [halves[0].split(',')[0], halves[0].split(',')[1], halves[1].split(',')[0], halves[1].split(',')[1]]
        try:
            taxa = [taxon_namespace[int(i) - 1] for i in [tokens[0], tokens[1], tokens[2], tokens[3]]]
        except Exception as e:
            print e
            print "Error reading line!"
            print line
            print tokens
            print len(taxon_namespace)
            exit()
        weight = 1.0
        return taxa, weight

    def run(self):
        dna = self.input_data['alignments'][0].clone(depth=2)

        taxon_namespace = dna.taxon_namespace

        del_list = []
        for t in taxon_namespace:
            if len([i for i in dna[t] if i != "-"]) == 0:
                del_list.append(t)

        for t in del_list:
            taxon_namespace.remove_taxon(t)
        
        print len(taxon_namespace)


        taxon_namespace.sort(key=lambda x: x.label)

        f = tempfile.NamedTemporaryFile(delete = False)

        dna.write_to_path(f.name + '.nex', 'nexus')
        f.flush()

        import subprocess

        t = time.time()

        pp = os.path.dirname(__file__) + '/parse_paup_qfile.sh'
        print "bash " + pp + " " + f.name
        out, err = subprocess.Popen(["bash", pp, f.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        print "Ran SVDQuartets in ", time.time() - t

        t = time.time()
        
        qfile = f.name + '.quartets'

        lines = open(qfile).readlines()

        quartets = WeightedQuartetSet(taxon_namespace)

        for l in lines:
            if len(l) and "discarded" not in l:
                (a,b,c,d), w = self.parse_line(l, taxon_namespace)
                quartets[(a,b,c,d)] = w

        print "read quartets in", time.time() - t

        self.result = {"quartets":quartets}
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
        if  tokens[-1] == "(rejected)":
            rejected = True
            tokens = tokens[:-1]
            
        if len(tokens) == 7:
            tokens = tokens[1:]

        try:
            taxa = [taxon_namespace[int(i) - 1] for i in [tokens[0], tokens[1], tokens[3], tokens[4]]]
        except:
            print "Error reading line!"
            print line
            print tokens
            print len(taxon_namespace)
            exit()
        weight = tokens[5]
        return taxa, weight


    def run(self):
        dna = self.input_data['alignments'][0].clone(depth=2)

        taxon_namespace = dna.taxon_namespace

        del_list = []
        for t in taxon_namespace:
            if len([i for i in dna[t] if i != "-"]) == 0:
                del_list.append(t)

        for t in del_list:
            taxon_namespace.remove_taxon(t)
        
        print len(taxon_namespace)


        taxon_namespace.sort(key=lambda x: x.label)

        f = tempfile.NamedTemporaryFile(delete = False)

        dna.write_to_path(f.name + '.nex', 'nexus')
        f.flush()

        import subprocess

        t = time.time()

        pp = os.path.dirname(__file__) + '/parse_paup.sh'
        print "bash " + pp + " " + f.name
        out, err = subprocess.Popen(["bash", pp, f.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        print "Ran SVDQuartets in ", time.time() - t

        t = time.time()
        

        lines = self.eliminate_header(out.split('\n'))
        lines = self.eliminate_footer(lines)
        
        print "cleaned up output in", time.time() - t

        t = time.time()

        quartets = WeightedQuartetSet(taxon_namespace)

        for l in lines:
            if len(l) and "discarded" not in l:
                (a,b,c,d), w = self.parse_line(l, taxon_namespace)
                quartets[(a,b,c,d)] = w

        print "read quartets in", time.time() - t

        self.result = {"quartets":quartets}
        return self.result


