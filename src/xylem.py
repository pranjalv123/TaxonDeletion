import argparse
import sys
import os
import importlib
import runpy
import Scheduler
import subprocess

template = """
#!/bin/bash
#PBS -l nodes=1:ppn=$ppn
#PBS -M pr@nj.al
#PBS -m a
#PBS -j oe
#PBS -l walltime=$walltime
#PBS -l naccesspolicy=singleuser
#PBS -q $queue
#PBS -t 0-$ntasks
#PBS -N $scriptname

export JOBID=`echo "$$PBS_JOBID" | cut -d"[" -f1`

export OMP_NUM_THREADS=10

echo $$PBS_JOBID
echo $$PBS_ARRAYID
echo $$JOBID

cd $$PBS_O_WORKDIR
mkdir -p  out.$$JOBID
cd  out.$$JOBID



for i in $$(seq $PBS_ARRAYID  $njobs $ntasks )
do
    python $scriptname  --nodot --single $$i --job $$JOBID
done
"""

class Xylem:
    def __init__(self):
        ap = argparse.ArgumentParser(description = "Run a Xylem script")
        ap.add_argument('command', type=str, help="Options: queue, single")
        ap.add_argument('scriptname', type=str, help="The script to run")

        ap.add_argument('--queue', '-q', type=str, help="The queue to use", dest='queue', default='secondary')
        ap.add_argument('--ppn', type=int)

        ap.add_argument('--njobs', type=int)
        ap.add_argument('--walltime', '-w', type=str, default='04:00:00', dest='walltime')
        

        self.args = vars(ap.parse_args())

        self.get_ntasks()
        qsub = self.gen_qsub()
        f = open(self.args['scriptname'] + '_qsub', 'w')
        f.write(qsub)

        subprocess.call(['qsub', f.name]
        
    def get_ntasks(self):
        sch = scheduler.SingleTaskScheduler(-1)
        mod = importlib.import(self.args.ap.script)
        mod.run(sch)
        self.args['ntasks'] = sch.total
        
    def gen_qsub(self):
        return string.Template(template).substitute(self.args)

    


def run_script():
    pass
