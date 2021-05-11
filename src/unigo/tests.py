from .utils import loadUniprotCollection
from . import Unigo as createGOTreeTest

DEFAULT_PROTEOME=f"{os.path.dirname(os.path.abspath(__file__))}/data/uniprot-proteome_UP000000807.xml.gz"
def run(nDummy=50, nTop=10, target='fisher'):
    print(f"Testing local implementation of {target}")
    
    print(f"Loading test proteome {DEFAULT_PROTEOME}")
    uTaxid, uColl  = loadUniprotCollection(DEFAULT_PROTEOME, strict=False)
    deltaUniprotID = generateDummySet(uColl, nDummy)
    
    print("Creatin unigo Tree")        
    unigoTree = createGOTreeTest(backgroundUniColl = uColl,
                                proteinList       = expUniprotID,
                                fetchLatest       = False)
    
    if target == 'fisher' or target == 'convert':
        print("Computing ORA")
        rankingsORA = unigoTree.computeORA(deltaUniprotID)
        print(f"Test Top - {nTop}\n{rankingsORA[:nTop]}")

    if arguments['convert']:
        print("Testing tree serialization")
        d = unigoTree.tree.f_serialize()
        print(d.asDict)