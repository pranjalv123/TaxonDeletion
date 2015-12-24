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
        if task.uid in self.scheduled:
            print "Didn't schedule", str(task), task.status()
            return
        #task.set_status("scheduled")
        self.scheduled.add(task.uid)
        print "Scheduling", task
        self.queue.append(task)
        self.tasks[task.uid] = task
        
    def add(self, plfun):
        self.pipelines.append(plfun(self))

    def run(self):
        while len(self.pipelines):
            p = self.pipelines.pop()
            self.current_pl = p
            p.ready()
            self.run_pl()
        
    def run_pl(self):
        print "Running scheduler"
        print "scheduled", self.scheduled
        print "running", self.running
        if '--memdebug' in sys.argv:
            gc.set_debug(gc.DEBUG_STATS | gc.DEBUG_INSTANCES | gc.DEBUG_OBJECTS)

        while len(self.queue):
            task = self.queue.popleft()

            gc.collect()


            try:
                task.execute(self.cache, self.regen)
                self.current_pl.todot()
                task.set_status("complete")
                for t2 in task.allows():
                    t2.req_complete(task)
                    print "allows", t2
                    if t2.status() == "ready":
                        print t2, "ready"
                        self.schedule(t2)
            except Task.DependenciesNotCompleteException:
                for dep in task.dependencies:
                    self.schedule(dep)
                print "Depends on:"
                for i in task.dependencies:
                    print str(i), i.status()
                self.queue.append(task)
                print
                print "QUEUE:"
                for i in self.queue:
                    print i
                print
            except Exception as e:
                open("ERRORS", 'a').write(str(task), 'failed with', e)
                print "FAILEDTASK", e
                self.queue.clear()
                break


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
                    
#if '--distributed' in sys.argv:
#    from mpi4py import MPI                    
                    
class MPIScheduler:
    def __init__(self, cache = True, regen=False, hostrank=0):
        print "ERROR: MPIScheduler is currently not working; try DistributedSerialScheduler"
        exit()
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.cache = cache
        self.regen = regen
        
        if self.rank != hostrank:
            JobRunner(self.comm, self.rank, self.cache, self.regen)
            exit()
        self.unready = []
        self.scheduled = []
        self.tasks = {}
        self.running = set()
        self.freepes = []

    def schedule(self, task):
        task.set_status("scheduled")
        self.scheduled.append(task)
        self.tasks[task.uid] = task
    def run(self):
        print "Running scheduler"
        print "scheduled", self.scheduled
        print "running", self.running
        print "freepes", self.freepes

        while len(self.scheduled) or len(self.running):
            pe, uid, result = self.comm.recv(source=MPI.ANY_SOURCE, tag=AVAIL)
            
            if uid:
                done_task = self.tasks[uid]
                done_task.set_status("complete")
                print "pe", pe, "finished", done_task
                            
                done_task.result = result
                self.running.remove(done_task)
                
                for t2 in done_task.allows():
                    t2.req_complete(done_task)
                    print "allows", t2
                    if t2.status() == "ready":
                        print t2, "ready"
                        self.schedule(t2)
            
            self.freepes.append(pe)
            toremove = set()
                    
            while len(self.freepes) and len(self.scheduled):
                task = self.scheduled.pop()
                if task.local:
                    print "Running", task, "locally"
                    
                    task.execute(self.cache, self.regen)
                    
                    task.set_status("complete")
                    for t2 in task.allows():
                        t2.req_complete(task)
                        print "allows", t2
                        if t2.status() == "ready":
                            print t2, "ready"
                            self.schedule(t2)

                else:
                    self.running.add(task)
                    self.execute(task, self.freepes.pop())
        print "Killing PEs"
        print self.freepes
        for pe in self.freepes:
            self.execute(Tasks.Exit(), pe)
    def execute(self, task, pe):
        print "sending", task, "to", pe
        self.comm.send(task, dest=pe)
            

class JobRunner:
    def __init__(self, comm, rank, cache, regen):
        print "Starting JobRunner", rank
        comm.send((rank, None, None), dest=0, tag=AVAIL)
        while True:
            task = comm.recv(source=0)
            print "Running", task.desc(), "on", rank
            if task.EXIT:
                return
            result = task.execute(cache, regen)
            comm.send((rank, task.uid, result), dest=0, tag=AVAIL)
