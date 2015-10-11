import Scheduler
import cPickle


class Pipeline:
    def __init__(self, scheduler):
        self.tasks = []
        self.scheduler = scheduler
    def add_task(self, task):
        self.tasks.append(task)
        task.pipeline = self
        return task
    def run(self, cache=True, regen=False):
        for task in self.tasks:
            if task.status() == "ready":
                self.scheduler.schedule(task)
        self.scheduler.run()
    
