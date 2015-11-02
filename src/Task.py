"""
Xylem - Phylogenetic Pipelines with MPI

Task.py contains routines for managing pipeline tasks.
    
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
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid
import time

class Task(object): #should be a "new style class" for inheritance purposes
    def __init__(self, *args, **kwargs):
        cachefile=None
        local=False
        cache=True
        if 'cachefile' in kwargs:
            cachefile = kwargs['cachefile']
        if 'local' in kwargs:
            local = kwargs['local']
        if 'cache' in kwargs:
            local = kwargs['cache']
        self.dependencies = set()
        self.depended = set()
        self.input_data = {}
        self._status="waiting"

        self.EXIT=False
        self.uid = uuid.uuid1()

        # flags that can be set by the user
        self.local = local
        self.cache = cache
        self.cachefile = cachefile

        self.setup(*args, **kwargs)

    #should be implemented here
    
    def require(self, otherTask): #don't run until otherTask has run
        self.dependencies.add(otherTask)
        otherTask.depended.add(self)
        return self
    
    def depends(self): #return a list of tasks this depends on
        return self.dependencies

    def allows(self): #return a list of tasks that depend on this
        return self.depended

    def req_complete(self, otherTask): #called when a task this depends on is complete

        print otherTask.outputs()
        task_outputs = set(otherTask.outputs())
        
        for item in self.inputs():
            name, tpe = item
            if type(tpe) == tuple:
                if name not in self.input_data:
                    self.input_data[name] = []
                if (name, tpe) in task_outputs:
                    self.input_data[name].extend(otherTask.get_results()[name])
                else:
                    for elem in tpe:
                        if (name, elem) in task_outputs:
                            self.input_data[name].append(otherTask.get_results()[name])
            else:
                if item in task_outputs:
                    self.input_data[name] = tpe(otherTask.get_results()[name])
        
        for i in self.dependencies:
            if i.status() != "complete":
                return
        self.set_status("ready")
        return

    def get_results(self):
        return self.result
                    
    def status(self): #return "waiting" if not ready, "ready", "scheduled", "running", or "complete"
        if len(self.dependencies) == 0 and self._status == "waiting":
            return "ready"
        else:
            return self._status

    def set_status(self, status):
        self._status = status

    def finish(self):
        self.set_status("complete")

    def set_inputs(self, inputs):
        self._inputs = inputs

    def storefile(self):
        if self.cachefile:
            return self.cachefile
        return ""
    
    def execute(self, cache=True, regen=False):
        t0 = time.clock()
        if not cache or not self.cache:
            self.result = self.run()
            print "Running", self, "took", time.clock() - t0, "seconds"
            return self.result
        filename = self.storefile()
        if not filename:
            print "Filename", filename, "not found", self.cache
            self.result =  self.run()
            print "Running", self, "took", time.clock() - t0, "seconds"
            return self.result
        try:
            if not regen:
                self.result = self.read(filename)
                if self.result:
                    print "Reading cache of", self, "took", time.clock() - t0, "seconds"
                    return self.result
        except Exception as e:
            print "Couldn't read file!", filename
            print e
        self.result = self.run()
        self.write(filename)
        if debug:
            print "Wrote", self, "to", filename
        for i in self.outputs():
            assert(i[0] in self.result.keys())
        print "Running", self, "took", time.clock() - t0, "seconds"
        return self.result

    def __hash__(self):
        try:
            return hash(self.uid)
        except:
            return super(Task, self).__hash__()
    def __eq__(self,other):
        try:
            return self.uid == other.uid
        except:
            return super(Task, self).__eq__(other)
    def id(self): 
         return self.uid

    
    #may be reimplemented by children
    def desc(self):
        return self.__repr__()
    def write(self, fname):
        f = open(fname, 'w')
        cPickle.dump(self.result, f, protocol=2)
    def read(self, fname):
        f = open(fname)
        return cPickle.load(f)
    def setup(self, *args, **kwargs):
        pass
    #need to be implemented by children
    def run(self): #actually do the thing
        pass
    def inputs(self): #returns a set of attributes it expects to see in its input
        pass
    def outputs(self): #returns a set of attributes it will have in its output
        pass
