import dendropy
import numpy as np
import sys
import argparse
from DataSet import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Delete taxa from PHYLIP sequence files")

#    parser.add_argument('-e', '--estimatedgenetrees', required=False)
    parser.add_argument('-t', '--truegenetrees', required=False)
    parser.add_argument('-a', '--alignments', required=True)
    parser.add_argument('-s', '--speciestree', required=True)
    parser.add_argument('-l', '--maxlen', type=int, required=True)
    parser.add_argument('-n', '--ndelete', type=int, required=True)
    parser.add_argument('-v', '--variance', type=int, required=True)
    parser.add_argument('-r', '--restrict', type=int, required=True)
    parser.add_argument('-o', '--output',  required=True)
    

    tn = dendropy.TaxonNamespace()
    
    args = vars(parser.parse_args())
    output = args['output']

    ndelete = args['ndelete']
    nrestrict = args['restrict']
    sigma = (args['variance'])**0.5
    if 'maxlen' in args:
        maxlen = args['maxlen']
    else:
        maxlen = None
        
    # if 'estimatedgenetrees' in args:
    #     estimatedgenetrees = dendropy.TreeList.get_from_path(args['estimatedgenetrees'], 'newick', taxon_namespace=tn)
    # else:
    #     estimatedgentrees = None
        
    if 'truegenetrees' in args and args['truegenetrees']:
        truegenetrees = dendropy.TreeList.get_from_path(args['truegenetrees'], 'newick', taxon_namespace = tn)
    else:
        truegenetrees = None
 
    print "Reading alignments..."
    alignments = read_multiphylip(args['alignments'], taxon_namespace = tn)

    print "Reading trees..."
    speciestree = dendropy.Tree.get_from_path(args['speciestree'], 'newick', taxon_namespace = tn)

    ds = DataSet(tn, truegenetrees, alignments, speciestree)

    print "Restricting taxa..."
    restricted = ds.delete_taxa_uniform(nrestrict, genetrees=True, seqs=True, speciestree=True, maxlen=maxlen)

    print "Deleting taxa..."    
    missing = restricted.delete_taxa(ndelete, sigma, genetrees=True, seqs=True, maxlen=maxlen)


    s = missing.seqs[0]
    
    estimated = DataSet(tn, None, missing.seqs, missing.speciestree)
    print "Estimating trees..."
    estimated.est_trees_fasttree()

    print "Writing data..."
    missing.write(output, speciestree=True)
    estimated.write(output + '_est', genetrees=True, seqs=True)

    
