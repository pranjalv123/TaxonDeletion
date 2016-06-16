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
import traceback

class DependenciesNotCompleteException(Exception):
    pass

class Task(object):
    cachefile_set = {}
    def __init__(self, *args, **kwargs):
        self.stackframe = traceback.extract_stack()[:-1]

        cachefile=kwargs.pop('cachefile', None)
        local=kwargs.pop('local', False)
        cache=kwargs.pop('cache', True)
        self.regen=kwargs.pop('regen', False)
        is_result=kwargs.pop('is_result', False)
                
        self.result = None
        self.dependencies = set()
        self.depended = set()
        self.input_data = {}
        self._status="waiting"

        self._dot_status="waiting"

        self.EXIT=False
        self.uid = uuid.uuid1()

        # flags that can be set by the user
        self.local = local
        self.cache = cache
        self.cachefile = cachefile
        self._is_result = is_result


        if cachefile and cache:
            if cachefile in Task.cachefile_set:
                print "ERROR", cachefile, "used as cache file for two tasks!:"
                print str(Task.cachefile_set[cachefile])
                print str(self)
                exit(-1)
        Task.cachefile_set[cachefile] = self

        
        
        self.setup(*args, **kwargs)


    #should be implemented here

    def is_result(self):
        if len(self.outputs()) == 0:
            return True
        return self._is_result

    #don't run until otherTask has run
    def require(self, *tasks):
        for otherTask in tasks:
            self.dependencies.add(otherTask)
            otherTask.depended.add(self)
        return self
    
    def depends(self): #return a list of tasks this depends on
        return self.dependencies

    def allows(self): #return a list of tasks that depend on this
        return self.depended

    def req_complete(self, otherTask): #called when a task this depends on is complete
        if self.status() == "complete":
            return
        print otherTask.outputs()
        task_outputs = set(otherTask.outputs())
        
        for item in self.inputs():
            if item == "*":
                for name, _ in task_outputs:
                    self.input_data[name] = otherTask.get_results()[name]
                continue
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
                    self.input_data[name] = tpe( otherTask.get_results()[name])                
            
        for i in self.dependencies:
            if i.status() != "complete":
                return
        self.set_status("ready")
        self._dot_status="ready"
        return

    def get_results(self):
        return self.result
                    
    def status(self): #return "waiting" if not ready, "ready", "scheduled", "running", or "complete"
        if len(self.dependencies) == 0 and self._status == "waiting":
            self._dot_status="ready"
            return "ready"
        else:
            return self._status

    def set_status(self, status):
        self._status = status

    def finish(self):
        self._dot_status="complete"
        self.set_status("complete")

    def set_inputs(self, inputs):
        self._inputs = inputs

    def storefile(self):
        if self.cachefile:
            return self.cachefile
        return ""

    def get_cache(self, cache=True, regen=False):
        cache &= self.cache
        regen |= self.regen
        if self.result:
            return True
        if cache and (not regen) and self.storefile():
            try:
                t0 = time.time()
                print "Trying to read cache of", self, "from", self.storefile()
                self.result = self.read(self.storefile())
                if self.result:
                    print "Reading cache of", self, "took", time.time() - t0, "seconds"
                    return True
            except Exception as e:
                print "Failed", e
                return False

    def write_cache(self, cache=True, regen=False):
        cache &= self.cache
        if cache and self.storefile():
            print "Writing to", self.storefile()
            try:
                self.write(self.storefile())
            except Exception, e:
                print "Error writing cache:"
                print e
        if self.storefile():
            print "should have written to", self.storefile()
        
    def execute(self, cache=True, regen=False):
        self.get_cache(cache, regen)
        if self.result or self.status() == "complete":
            return self.result

        if not self.status() == "ready":
            raise DependenciesNotCompleteException
        
        print "Starting", str(self)
        self._dot_status="running"
        t0 = time.time()
        self.result = self.run()
        print "Running", self, "took", time.time() - t0, "seconds"
        
        self.write_cache()
    
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
    def __str__(self):
        
        return '(' + self.__class__.__name__ + '; ' + str(self.desc()) + ')'

    def desc(self):
        return str(self.uid)
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
        return {}
    def inputs(self): #returns a set of attributes it expects to see in its input
        return []
    def outputs(self): #returns a set of attributes it will have in its output
        return []
