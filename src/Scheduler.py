"""
Xylem - Phylogenetic Pipelines with MPI

Scheduler.py contains routines for distributing tasks across worker 
nodes with MPI.
    
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

import Task
#import Tasks
import sys
import gc
import collections


class Exit(Task.Task):
    def setup(self, *args, **kwargs):
        self.EXIT = True
    def desc(self):
        return ""
    def inputs(self):
        return []
    def outputs(self):
        return []


def Scheduler(cache=True, regen=False):
    if '--nocache' in sys.argv:
        cache=False
    if '--regen' in sys.argv:
        regen=True    
    if '--serial' in sys.argv:
        return SerialScheduler(cache, regen)
    if '--distributed' in sys.argv:
        return DistributedSerialScheduler(cache, regen)
    if '--single' in sys.argv:
        return SingleTaskScheduler(int(sys.argv[sys.argv.index('--single') + 1]), cache, regen)
    return MPIScheduler()

AVAIL=1

class SerialScheduler:
    def __init__(self, cache=True, regen=False):
        self.cache = cache
        self.regen = regen
        self.unready = []
        self.queue = collections.deque()
        self.tasks = {}
        self.running = set()
        self.pipelines = []
        self.scheduled = set()
    def schedule(self, task):
        print task
        print task.uid
        print task.status()
        if task.uid in self.scheduled:
            return False
        #task.set_status("scheduled")
        self.queue.append(task)
        self.tasks[task.uid] = task
        self.scheduled.add(task.uid)
        return True
        
    def add(self, plfun):
        self.pipelines.append(plfun(self))

    def run(self):
        while len(self.pipelines):
            p = self.pipelines.pop()
            self.current_pl = p
            p.ready()
            self.run_pl()
        
    def run_pl(self):

        if '--memdebug' in sys.argv:
            gc.set_debug(gc.DEBUG_STATS | gc.DEBUG_INSTANCES | gc.DEBUG_OBJECTS)

        while len(self.queue):
            task = self.queue.popleft()
            try:
                print "trying", task
                
                task.execute(self.cache, self.regen)
                task.set_status("complete")
                for t2 in task.allows():
                    t2.req_complete(task)
                    print "allows", t2
                    if t2.status() == "ready":
                        print t2, "ready"
#                        if self.schedule(t2):
#                            print "because we finished", str(task)
                self.current_pl.todot()
                
            except Task.DependenciesNotCompleteException:
                for dep in task.dependencies:
                    if self.schedule(dep):
                        print "to enable", str(task)
                print "ADDING", task.uid, "TO QUEUE"
                self.queue.append(task)


import os
import sys
def rank_hack():
    if 'ALPS_APP_PE' in os.environ:
        return int(os.environ['ALPS_APP_PE']), int(os.environ['NUM_PES'])
    if 'OMPI_COMM_WORLD_RANK' in os.environ:
        return int(os.environ['OMPI_COMM_WORLD_RANK']), int(os.environ['OMPI_COMM_WORLD_SIZE'])
    return 0, 1
                
class DistributedSerialScheduler:
    def __init__(self, cache = True, regen = False):
        self.sched = SerialScheduler(cache, regen)
        self.index = 0
#        self.comm = MPI.COMM_WORLD
        self.rank,self.size = rank_hack()#self.comm.Get_rank()

        sys.stdout = open(str(self.rank)  + '_out', 'w')
        sys.stderr = open(str(self.rank)  + '_err', 'w')

        print "RANK,SIZE:", self.rank,self.size
#        self.size = self.comm.Get_size()
        self.pipelines = []
    def add(self, plfun):
        if self.index == self.rank:
            self.pipelines.append(plfun(self.sched))
        self.index += 1
        self.index %= self.size
    def run(self):
        for p in self.pipelines:
            self.sched.current_pl = p
            p.ready()
            self.sched.run_pl()


class SingleTaskScheduler:
    def __init__(self, n, cache = True, regen = False):
        self.sched = SerialScheduler(cache, regen)
        self.index = n
        self.pipeline = None
        self.total = 0
    def add(self, plfun):
        if self.index == self.total:
            self.pipeline = (plfun(self.sched))
        self.total += 1
    def run(self):
        print "JOB:", self.index, "/", self.total
        self.sched.current_pl = self.pipeline
        self.pipeline.ready()
        self.sched.run_pl()


