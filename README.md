# Xylem: Phylogenetic Pipelines with MPI

Xylem is a library for easily constructing and running phylogenetic
pipelines. The initial scope of Xylem is limited to method evaluation,
but it could also be used for analyzing biological datasets.

Xylem works on any MPI-enabled system, and automatically schedules and
distributes jobs across nodes.

# Tasks

The building block of a Xylem pipeline is a Task. A Task represents
a transformation of input data to output data based on some
parameters.


## Defining a new task

New Tasks are defined by subclassing the Task class.

A Task subclass *MUST* define three functions: inputs(), outputs(),
and run().

1. inputs() returns a (possibly empty) list of tuples, each with type
(str, type). These represent the input data for the task. The string
is the name of the input, and the type is the type of the input. If
the input should be a collection of items, don't use list as a type;
rather, give a tuple of types, where each item can be of any of those
types. 

2. output() also returns a (possibly empty) list of tuples, each with
type (str, type). These represent the output data for the task.

3. run() does the task. The task runs when all the input data is ready,
and is provided in the dictionary self.input_data, which can be
indexed by the names given in inputs().

You may also override the write(self, fname) and read(self, fname)
functions. These should, respectively, write and read the self.result
dictionary and are used for caching results. By default these pickle
the dictionary.

Furthermore, you may override the setup function to provide non-input
parameters for tasks.

There are also three flags you can set in the constructor or when you
instantiate the object:

1. local (default: False): If true, the task will be run on the master
processor, and will not be distributed to a worker with MPI. This is
useful for quick tasks.

2. cache (default: True): If true, the output will be cached; if the
cache file exists, the task will not be run and will only be read from
the cache file

3. cachefile (default: ""): a path to a file used for the
cache. Caching only happens if cache=True and cachefile is a non-empty
string.

This is useful not only for the standard purposes of not unnecessarily
re-running jobs, but also for running jobs on systems with limited
wallclock times, like Blue Waters, effectively acting as a
checkpointing system.

# Pipeline

A pipeline is a group of tasks linked by dependencies in a directed
acyclic graph (DAG).

Tasks can be added to a pipeline with Pipeline.add_task(task).

Dependencies can be added to a task using Task.require(otherTask).

Both these functions have the feature of returning the added or
depending task, so if you have a pipeline pl, you can create tasks,
add them to the pipeline, and add dependencies in just one line

	task1 = pl.add_task(MyTask(...args...))
	task2 = pl.add_task(MyOtherTask(...args...)).require(task1)

or if you want to be really fancy, describe an entire pipeline or a
section of a pipeline in a single line:

	pl.add_task(MyOtherTask(...args...)).require(pl.add_task(MyTask(...args...)))

but this can get hard to read so is not necessarily recommended.

# Scheduler

A scheduler runs tasks with MPI. Users typically don't need to
interact with the scheduler, but you can have multiple pipelines with
the same scheduler.

When creating a scheduler, you can turn off caching globally with by
setting cache=False in the constructor. You can also force disable
reading from the cache, but allow new cached results to be written by
setting regen=True.

# Future improvements

Currently, the software does not check for inconsistencies in the
pipeline. For example, if there are namespace collisions, there will
be race conditions; if some task does not have every input as the
outputs of one of the dependencies, the pipeline will never finish; if
the dependency graph has cycles, the pipeline may never finish.

Sometimes, outputs of cached jobs will be loaded even if they're not
necessary. 