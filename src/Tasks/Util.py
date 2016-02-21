"""
Xylem - Phylogenetic Pipelines with MPI

Util.py contains useful functions that handle minor tasks but
typically don't affect data much.
    
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

import dendropy
import xylem.Task

class CastName(xylem.Task):
    def setup(self, inname, outname, tpe, *args, **kwargs):
        self.local=True
        self.cache=False
        self.tpe = tpe
        self.inname = inname
        self.outname = outname
    def desc(self):
        return self.inname + ' to ' + self.outname
    def inputs(self):
        return [(self.inname, self.tpe)]
    def outputs(self):
        return [(self.outname, self.tpe)]
    def run(self):
        self.result =  {self.outname : self.input_data[self.inname]}
        return self.result        

class Head(xylem.Task):
    def setup(self, n, *args, **kwargs):
        self.n = n
    def desc(self):
        return str(self.n)
    def inputs(self):
        return [('genetrees', dendropy.TreeList)]
    def outputs(self):
        return [('genetrees', dendropy.TreeList)]
    def run(self):
        gt = self.input_data["genetrees"]
        self.result = {'genetrees': gt[:min(self.n, len(gt))]}
        return self.result

class GtreesToStree(xylem.Task):
    def inputs(self):
        return [('genetrees', dendropy.TreeList)]
    def outputs(self):
        return [('estimatedspeciestree', dendropy.Tree)]
    def run(self):
        self.result = {'estimatedspeciestree':self.input_data['genetrees'][0]}
        return self.result
    
    
class Append(xylem.Task):
    def inputs(self):
        return [('genetrees', dendropy.TreeList), ('genetrees_2', dendropy.TreeList)]
    def outputs(self):
        return [('genetrees', dendropy.TreeList)]
    def run(self):
        gt = self.input_data["genetrees"]
        for t in self.input_data["genetrees_2"]:
            gt.append(t)
        self.result = {'genetrees': gt}
        return self.result


class RemoveBigPolytomies(xylem.Task):
    def setup(self, n):
        self.n = n
    def inputs(self):
        return [('genetrees', dendropy.TreeList)]
    def outputs(self):
        return [('genetrees', dendropy.TreeList)]
    def run(self):
        tl = self.input_data['genetrees']
        output = dendropy.TreeList(taxon_namespace = tl.taxon_namespace)
        for t in tl:
            if max([len(i.child_nodes()) for i in t.internal_nodes()]) < self.n:
                output.append(t)
        self.result = {'genetrees':output}
        return self.result

class ResolveBigPolytomies(xylem.Task):
    def setup(self, n):
        self.n = n
    def inputs(self):
        return [('genetrees', dendropy.TreeList)]
    def outputs(self):
        return [('genetrees', dendropy.TreeList)]


    def run(self):
        tl = self.input_data['genetrees']
        output = dendropy.TreeList(taxon_namespace = tl.taxon_namespace)
        
        for t in tl:
            t.resolve_polytomies(limit=self.n)
        self.result = {'genetrees':tl}
        return self.result


class Const(xylem.Task):
    def setup(self, name, val, tpe=None):
        self.val = val
        self.tpe = type(val)
        self.name = name
        if tpe:
            self.tpe = tpe
    def inputs(self):
        return []
    def outputs(self):
        return [(self.name, self.tpe)]
    def run(self):
        self.result = {self.name:self.val}
        return self.result

    

class Echo(xylem.Task):
    def setup(self, desc):
        self.describe = desc
    def inputs(self):
        return ['*']
    def run(self):
        print self.describe
        for name in self.input_data:
            print name, '=', self.input_data[name]
            print
