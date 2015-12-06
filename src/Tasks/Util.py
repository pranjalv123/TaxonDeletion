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
    
    
    
