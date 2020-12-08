from .tree import setOntology, createGoTree
from . import stat_utils

"""
Minimal ontology tree of a uniprot collection

Merging Uniprot Object Collection 
Consrtuctor

"""

def prune(unigoObj, predicate):
    """Create a Unigo object by Pruning the supplied one

        Parameters
        ----------
        unigo : 
        predicate:
    """
    _tree = unigoObj.drop(predicate)
    return Unigo(previous=_tree)

class Unigo:
    def __init__(self, previous=None, backgroundUniColl=None, proteinList=None, owlFile=None, ns="biological process"):
        """Create a Unigo object

        Parameters
        ----------
        owlFile : path to the ontology owl file.
        backgroundUniColl : collection of uniprot objects, defining the background population
        proteinList : A list of uniprot identifiers
        [Options]
        ns : a subset of ontology, default:biological process

        TODO: Check if all protein lsit are members of bkgUniCol

        """
        if previous:
            self.tree = previous
            return
        
        if backgroundUniColl is None or  proteinList is None:
            raise ValueError(f"parameters backgroundUniColl and proteinList are required")
        try :
            if owlFile is None:
                print("Fetching ontology")
                setOntology(url="http://purl.obolibrary.org/obo/go.owl")
            else:
                setOntology(file=owlFile)
        except Exception as e:
            print(f"Could not create ontology")
            print(e)
        
        self.tree = createGoTree(ns=ns, proteinList=proteinList, uniprotCollection=backgroundUniColl)

# Define rich view of stat results for notebook ?
def dumpStat():
    pass

#import json

#res = {}

#for node in goTreeObjExp.walk():
#    if node.pvalue:
#         res[ node.name ] = {
#             "name"          : node.name,
#             "pvalue"        : node.pvalue,
#             "proteineTotal" : node.getMembers(nr=True),
#             "proteineSA"    : list ( set (saList) & set(node.getMembers(nr=True)) )
#         }
        
#with open('TP_ORA.json', 'w') as fp:
#    json.dump(res, fp)
#res