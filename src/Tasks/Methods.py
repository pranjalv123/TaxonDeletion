"""
Xylem - Phylogenetic Pipelines with MPI

Methods.py contains routines for interacting with phylogenetic methods.
    
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

import xylem.Task
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid
import ASTRID
import dendropy
import tempfile
import os

class RunFastTree(xylem.Task):
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("genetrees", dendropy.TreeList)]
    
    def desc(self):
        return ""
    
    def run(self):
        self.seqs = self.input_data["alignments"]
        sio = StringIO.StringIO()
        for seq in self.seqs:
            seq.write_to_stream(sio, schema="phylip", suppress_missing_taxa=True)
        proc = subprocess.Popen(['fasttree', '-nt', '-gtr', '-nopr', '-gamma', '-n', str(len(self.seqs))], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        trees, err = proc.communicate(sio.getvalue())
        print err
        genetrees = dendropy.TreeList.get_from_string(trees, 'newick')
        self.result = {"genetrees":genetrees}
        return self.result
    
class RunASTRID(xylem.Task):
    def setup(self, distmethod=None, *args, **kwargs):
        self.distmethod = distmethod
        if self.distmethod==None:
            self.distmethod="auto"
    def desc(self):
        return self.distmethod
    
    def inputs(self):
        return [("genetrees", dendropy.TreeList)]
    def outputs(self):
        return [("estimatedspeciestree", dendropy.Tree)]
    def run(self):
        a = ASTRID.ASTRID(self.input_data["genetrees"])
        print a
        print self.input_data["genetrees"]
        a.run(self.distmethod)
        self.result = {"estimatedspeciestree": a.tree}
        print a.tree
        print type(a.tree)
        del(a)
        return self.result

class RunASTRAL(xylem.Task):
    def setup(self, *args, **kwargs):
        pass
    def desc(self):
        return ""
    def inputs(self):
        return [("genetrees", dendropy.TreeList)]
    def outputs(self):
        return [("estimatedspeciestree", dendropy.Tree)]
    def run(self):
        print "RUNNING ASTRAL"
        f = tempfile.NamedTemporaryFile()
        gt = self.input_data["genetrees"]
        gt = dendropy.TreeList([i for i in gt if len(i.leaf_nodes()) > 3])
        gt.write(path=f.name, schema='newick')
        print "ASTRAL", "-i", f.name
        proc = subprocess.Popen(['ASTRAL', '-i', f.name], stdout=subprocess.PIPE)
        streestr, err = proc.communicate()
        print err
        print streestr
        del(proc)
        stree = dendropy.Tree.get_from_string(streestr, 'newick')
        self.result = {"estimatedspeciestree": stree}
        return self.result

import Quartets
    
class RunWQMC(xylem.Task):
    def inputs(self):
        return [("quartets", Quartets.WeightedQuartetSet)]
    def outputs(self):
        return [("estimatedspeciestree", dendropy.Tree)]

    def run(self):
        f = tempfile.NamedTemporaryFile(delete=False)
        o = tempfile.NamedTemporaryFile(delete=False)
        ix = self.input_data["quartets"].write(f, 'wqmc', translate=True)
        print f.name
        f.close()
        proc = subprocess.Popen(['wQMC', '-weight', 'qrtt='+f.name, 'otre='+o.name])
        proc.wait()
        stree = dendropy.Tree.get_from_path(o.name, 'newick')
        for t in stree.taxon_namespace:
            t.label = ix[int(t.label)].label
        self.result = {"estimatedspeciestree":stree}
        return self.result


class RunWastral(xylem.Task):
    def setup(self, criterion='dp', score=False, exact=False, maximize=True):
        self.criterion = {'dp':'DPTripartitionScorer', 'bs':'BryantSteelTripartitionScorer',
                          'rf':'RFTripartitionScorer'}[criterion]
        self.score = score
        self.exact = exact
        self.maximize = maximize
        self.inputs_ = set()
        self.outputs_ = set()

        if criterion == 'dp' or criterion == 'bs':
            self.inputs_.add(("quartets", Quartets.WeightedQuartetSet))
        elif criterion == 'rf':
            self.inputs_.add(("genetrees", dendropy.TreeList))

        if not (exact or score):
            self.inputs_.add(("genetrees", dendropy.TreeList))
        if score:
            self.inputs_.add(("estimatedspeciestree", dendropy.Tree))
            self.outputs_.add(("score", float))
        else:
            self.outputs_.add(("estimatedspeciestree", dendropy.Tree))
    def desc(self):
        d = self.criterion
        if self.score:
            d += ' score '
        if self.exact:
            d += ' exact '
        if self.maximize:
            d += ' max '
        return d
    def inputs(self):
        return list(self.inputs_)
    def outputs(self):
        return list(self.outputs_)

    def run(self):

        o = tempfile.NamedTemporaryFile()
        
        args = ['wASTRAL', '-c', self.criterion, '-a', '~/.local/lib/astral.4.7.8.jar',  '-o', o.name]
        if self.maximize:
            args += ['--maximize']
        else:
            args += ['--minimize']

        if ("quartets", Quartets.WeightedQuartetSet) in self.inputs():
            qf = tempfile.NamedTemporaryFile(delete=False )
            self.input_data["quartets"].write(qf, 'wqmc')
            args += ['-q', qf.name]
            qf.flush()
        if ("genetrees", dendropy.TreeList) in self.inputs():
            gf = tempfile.NamedTemporaryFile(delete=False )
            self.input_data["genetrees"].write_to_path(gf.name, 'newick', suppress_edge_lengths=True)
            args += ['-g', gf.name]
            gf.flush()
        elif ("quartets", Quartets.WeightedQuartetSet) in self.inputs():
            gf = tempfile.NamedTemporaryFile(delete=False)
            gt = dendropy.simulate.star_tree(self.input_data["quartets"].tn)
            gt.resolve_polytomies()
            gt.write_to_path(gf.name, 'newick', suppress_edge_lengths=True)
            args += ['-g', gf.name]
            gf.flush()
        if self.score:
            gf = tempfile.NamedTemporaryFile(delete=False)
            print self.inputs()
            print self.input_data
            gt = self.input_data["estimatedspeciestree"]
            gt.write_to_path(gf.name, 'newick', suppress_edge_lengths=True)
            print "SCORING TREE:", str(gt)
            print (open(gf.name).read())
            args += ['-s', gf.name]
            
            gf.flush()
        
        print ' '.join(args)
        print args
        proc = subprocess.Popen(args)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.wait()
#        out, err = proc.communicate()
#        print out
#        print err

        if not self.score:
            stree = dendropy.Tree.get_from_path(o.name, 'newick')
            self.result = {"estimatedspeciestree":stree}
        else:
            score = float(open(o.name, 'r').read())
            self.result = {"score":score}
        return self.result

class RunSVDQuartets(xylem.Task):
    def setup(self):
        pass
    def inputs(self):
        return [("alignments", (dendropy.DnaCharacterMatrix,))]
    def outputs(self):
        return [("estimatedspeciestree", dendropy.Tree)]
    
    def run(self):
        
        dna = self.input_data['alignments'][0].clone(depth=2)

        
        taxon_namespace = dna.taxon_namespace

        del_list = []
        for t in taxon_namespace:
            if len(dna[t]) == 0:
                del_list.append(t)

        for t in del_list:
            taxon_namespace.remove_taxon(t)
        
        taxon_namespace.sort(key=lambda x: x.label)


        print "SVDQuartets"
        print len(dna)
        print len(dna[0])
        f = tempfile.NamedTemporaryFile()
        print f.name
        dna.write_to_path(f.name + '.nex', 'nexus')

        import subprocess

        pp = os.path.dirname(__file__) + '/parse_paup_svdtree.sh'
        print " ".join(["bash", pp, f.name])
        out, err = subprocess.Popen(["bash", pp, f.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        
        try:
            tree = dendropy.Tree.get_from_path(f.name + '.svdtree','newick')
        except IOError as e:
            print e
            print out, err
            exit(0)
        
        self.result = {"estimatedspeciestree":tree}
        return self.result


class RunMulRF(xylem.Task):
    def setup(self, niter=10):
        self.niter = niter
    def inputs(self):
        return [("genetrees", dendropy.TreeList)]
    def outputs(self):
        return [("estimatedspeciestree", dendropy.Tree)]    
    def run(self):
        f = tempfile.NamedTemporaryFile()
        gt = self.input_data["genetrees"]
        gt = dendropy.TreeList([i for i in gt if len(i.leaf_nodes()) > 3])
        gt.write(path=f.name, schema='newick')
        bestscore = 9999999999999
        stree = None
        
        for i in range(self.niter):
            st = tempfile.NamedTemporaryFile()
            args = ['MulRFSupertree', '-i', f.name, '-o', st.name]
            proc = subprocess.Popen(args)
            proc.wait()
            lines = st.readlines()
            print "LINES:"
            for l in lines:
                print l
            score = float(lines[0][len("[ Species Tree: Unrooted RF Score = "):-2])
            tree = lines[1]
            print score
            if score < bestscore:
                bestscore = score
                stree = dendropy.Tree.get_from_string(lines[1], 'newick')
            
            
        self.result = {"estimatedspeciestree":stree}
        return self.result
