import Task
import DataSet
import cPickle
import numpy as np
import subprocess
import StringIO
import uuid

class CompareTrees(Task):
    def __init__(self, *args, **kwargs):
        super(CompareTrees, self).__init__(*args, **kwargs)
    def inputs(self):
        return[("truespeciestree", dendropy.Tree), ("estimatedspeciestree", dendropy.Tree)]
    def outputs(self):
        return[("rfdistance", float)]
    def run(self):
        truetree = self.input_data["truespeciestree"]
        estimatedtree = self.input_data["estimatedspeciestree"]
        tn = estimatedtree.taxon_namespace
        truetree.migrate_taxon_namespace(tn)
        diff = dendropy.calculate.treecompare.false_positives_and_negatives(truetree, estimatedtree)[0]
        self.result ={"rfdistance":diff}
        print "RF distance is", diff
        return self.result
