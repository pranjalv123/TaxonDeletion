import dendropy
import sys
dendropy.datamodel.charmatrixmodel.DnaCharacterDataSequence = dendropy.DnaCharacterMatrix.DnaCharacterDataSequence
dendropy.datamodel.charmatrixmodel.RnaCharacterDataSequence = dendropy.RnaCharacterMatrix.RnaCharacterDataSequence
dendropy.datamodel.charmatrixmodel.NucleotideCharacterDataSequence = dendropy.NucleotideCharacterMatrix.NucleotideCharacterDataSequence
dendropy.datamodel.charmatrixmodel.ProteinCharacterDataSequence = dendropy.ProteinCharacterMatrix.ProteinCharacterDataSequence
dendropy.datamodel.charmatrixmodel.RestrictionSitesCharacterDataSequence = dendropy.RestrictionSitesCharacterMatrix.RestrictionSitesCharacterDataSequence
dendropy.datamodel.charmatrixmodel.InfiniteSitesCharacterDataSequence = dendropy.InfiniteSitesCharacterMatrix.InfiniteSitesCharacterDataSequence
dendropy.datamodel.charmatrixmodel.StandardCharacterDataSequence = dendropy.StandardCharacterMatrix.StandardCharacterDataSequence

from Methods import *
from Delete import *
from Phylo import *
from Readers import *
from Util import *
from Writers import *
from Quartets import *
from Sequence import *
from Bootstrap import *


