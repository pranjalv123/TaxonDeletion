"""
Xylem - Phylogenetic Pipelines with MPI

Pipeline.py contains routines for managing dependency graphs.
    
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

import Scheduler
import cPickle
import sys
from collections import defaultdict
import subprocess
import traceback

colormap = defaultdict(lambda : "black")
colormap['waiting']='red'
colormap['ready']='yellow'
colormap['running']='green'
colormap['complete']='blue'


import os
def stylemap(task):
    if not task.cache or not task.cachefile:
        return 'dashed'
    if os.path.isfile(task.cachefile):
        return 'solid'
    return 'dotted'


class Pipeline:
    def __init__(self, scheduler, desc):
        self.tasks = []
        self.task_set = set()
        self.scheduler = scheduler
        self.desc = desc
    def add_task(self, task):
        if task in self.task_set:
            return task

        self.tasks.append(task)
        self.task_set.add(task)
        return task
    def verify(self):
        for t in self.tasks:
            if t.outputs() == None:
                print t, "has invalid outputs()"
                print traceback.format_list(t.stackframe)
            if t.inputs() == None:
                print t, "has invalid inputs()"
                print traceback.format_list(t.stackframe)
            for i in t.inputs():
                good = False
                for d in t.depends():
                    if i in d.outputs() or i == '*':
                        good = True
                if not good:
                    print "Dependency not found:", t, i
                    print '\n'.join([':'.join([str(d), str(d.outputs())]) for d in t.depends()])
                    print ''.join(traceback.format_list(t.stackframe))
                    return False
        return True
    def prune(self):
        pass
    def todot(self):
        fname = self.desc + '.dot'
        f = open(fname, 'w')
        f.write('digraph pipeline {\n')
        for task in self.tasks:
            f.write('n'+task.uid.hex + '[label="' + str(task) + '",color="' + colormap[task._dot_status] + '",style="' + stylemap(task) + '"];\n')
            for dep in task.dependencies:
                f.write('n'+dep.uid.hex + '->' + 'n'+task.uid.hex +';\n')
        f.write('}\n')
        f.close()
        if '--nodot' not in sys.argv:
            subprocess.Popen(['dot', fname, '-Tpng', '-O'])
    def ready(self, cache=True, regen=False):
        self.todot()
        if not self.verify():
            return False
        for task in self.tasks:
            if task.is_result():
                self.scheduler.schedule(task)
