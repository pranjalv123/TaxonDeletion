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
import Tasks
from mpi4py import MPI

AVAIL=1

class Scheduler:
    def __init__(self, cache = True, regen=False):
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.cache = cache
        self.regen = regen
        
        if self.rank > 0:
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
            print
            print
            

            
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
                    if self.cache:
                        task.run_cached(self.regen)
                    else:
                        task.run()
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
            if cache:
                result = task.run_cached(regen)
            else:
                result = task.run()
            comm.send((rank, task.uid, result), dest=0, tag=AVAIL)
