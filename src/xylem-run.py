import argparse
import sys
import os
import imp
import runpy
import Scheduler
import subprocess32

import string

template = """#!/bin/bash
#PBS -l nodes=1:ppn=$ppn
#PBS -M pr@nj.al
#PBS -m a
#PBS -j oe
#PBS -l walltime=$walltime
#PBS -l naccesspolicy=singleuser
#PBS -q $queue
#PBS -N $scriptname
#PBS -t 0-$njobs
#PBS -o logs


export JOBID=`echo "$$PBS_JOBID" | cut -d"[" -f1`

export OMP_NUM_THREADS=1

echo $$PBS_JOBID
echo $$PBS_ARRAYID
echo $$JOBID

cd $$PBS_O_WORKDIR
mkdir -p  out.$$JOBID
mkdir -p  out.$$JOBID/logs/
cd  out.$$JOBID

cp ../$scriptname .

for i in $$(seq $$PBS_ARRAYID  $njobs $ntasks )
do
    python ../$scriptname  --nodot --single $$i --job $$JOBID $extra_args
done

mv $$LOGNAME logs/

"""

class Xylem:
    def __init__(self):
        ap = argparse.ArgumentParser(description = "Run a Xylem script")

        ap.add_argument('scriptname', type=str, help="The script to run")

        ap.add_argument('--queue', '-q', type=str, help="The queue to use", dest='queue', default='secondary')
        ap.add_argument('--ppn', type=int, default=1)

        ap.add_argument('--njobs', type=int, default=100)
        ap.add_argument('--walltime', '-w', type=str, default='04:00:00', dest='walltime')
        

        args, otherargs = ap.parse_known_args()

        self.args = vars(args)
        self.args['extra_args'] = ' '.join(otherargs)

        self.get_ntasks()
        qsub = self.gen_qsub()
        f = open(self.args['scriptname'] + '_qsub', 'w')
        print f.name
        f.write(qsub)
        args = 'qsub ' + f.name

        print "RUNNING", args
        os.system(args)

        os.system('env')
        
        
    def get_ntasks(self):
        sch = Scheduler.SingleTaskScheduler(-1)
        mod = imp.load_module("pl", open(self.args['scriptname']), './' + self.args['scriptname'], ('.py', 'U', 1))
        mod.run(sch)
        self.args['ntasks'] = sch.total
        self.args['njobs'] = min(self.args['njobs'], self.args['ntasks'])
        
    def gen_qsub(self):
        return string.Template(template).substitute(self.args)

    

if __name__ == "__main__":
    Xylem()
