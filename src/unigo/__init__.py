from .tree import setOntology, createGoTree, load
from . import stat_utils
import os
import json

"""
Minimal ontology tree of a uniprot collection

Merging Uniprot Object Collection 
Consrtuctor

"""

def prune(unigoObj, predicate):
    """Create a Unigo object by Pruning the supplied Unigo instance

        Parameters
        ----------
        unigoObj : 
        predicate: predicate function applied to prodived Unigo to drop or keep its elements
    """
    _tree          = unigoObj.tree.drop(predicate)
    _tree_universe = unigoObj.tree_universe.drop(predicate)
    return Unigo( previous=(_tree, _tree_universe) )


def chop(unigoObj, name=None, ID=None):
    """Create a Unigo object with its root extracted from the provided Unigo instance

        Parameters
        ----------
        unigoObj : 
        name: name of the node the look for in provided unigoObj. It will be the root of the returned Unigo
        ID: ID of the node the look for in provided unigoObj. It will be the root of the returned Unigo
    """
    
    _tree          = unigoObj.tree.newRoot(name=name, ID=name)
    _tree_universe = unigoObj.tree_universe.newRoot(name=name, ID=name)
    return Unigo( previous=(_tree, _tree_universe) )

def vloads(data):
    """Deserialize tree as dict structure fetched from api
    """
    _data = json.loads(data)
    return  Univgo(serial=_data)
 
    # Update leaf count
class Univgo:
    def __init__(self, serial = None, 
                       uniColl = None, 
                       owlFile=None, fetchLatest=True,
                       ns="biological process"):#, **kwargs):
        """Create a Univgo object, made to manipulate whole proteome ontology

        Parameters
        ----------
        owlFile : path to the ontology owl file.
        backgroundUniColl : collection of uniprot objects, defining the background population
        [Options]
        ns : a subset of ontology, default:biological process

        """
        if serial:
            self.single_tree = load(serial)
            return
        
        if uniColl is None :
            raise ValueError(f"parameters uniColl required")
        try :
            if owlFile is None and fetchLatest:
                print("Fetching ontology")
                setOntology(url="http://purl.obolibrary.org/obo/go.owl")#, **kwargs)
            else:
                print("Using package release ontology")
                setOntology(owlFile=f"{os.path.dirname(os.path.abspath(__file__))}/data/go.owl")
        except Exception as e:
            print(f"Could not create ontology")
            print(e)
        
        self.single_tree = createGoTree(         ns = ns,
                                  proteinList       = [ k for k in uniColl.keys() ], 
                                  uniprotCollection = uniColl)

    def serialize(self):
        return self.single_tree.f_serialize().asDict
    @property
    def dimensions(self):
        return self.single_tree.dimensions

def uloads(data, proteinIDList):
    """Generate a Unigo tree Object
    """
    _data=json.loads(data)

    tree_universe = Univgo(serial=_data).single_tree
    # Set predicate function to retain node in universal tree whih are path to provided proteins
    def keep_proteinID(node):
        return set(node.getMembers()) & set(proteinIDList)
    
    # Obtaine the tree
    exp_tree = tree_universe.drop(keep_proteinID, noCollapse=True, noLeafCountUpdate=True)
    # Reset eTags to keep only provided proteinIDList
    for n in exp_tree.walk():
        n.eTag = list(  set(proteinIDList) & set(n.eTag) )
    # update tree topology and leaf counts
    tree = exp_tree.collapse()
    exp_tree.leafCountUpdate()

    return Unigo(previous=(tree, tree_universe))

class Unigo:
    def __init__(self, previous = None, 
                       backgroundUniColl = None, 
                       proteinList       = None, 
                       owlFile=None, fetchLatest=True,
                       ns="biological process"):#, **kwargs):
        """Create a Unigo object

        Parameters
        ----------
        owlFile : path to the ontology owl file.
        backgroundUniColl : collection of uniprot objects, defining the background population
        proteinList : A list of uniprot identifiers
        [Options]
        ns : a subset of ontology, default:biological process

        TODO: Check if all protein list are members of bkgUniCol

        """
        if previous:
            self.tree, self.tree_universe = previous
            return
        
        if backgroundUniColl is None or  proteinList is None:
            raise ValueError(f"parameters backgroundUniColl and proteinList are required")
        try :
            if owlFile is None and fetchLatest:
                print("Fetching ontology")
                setOntology(url="http://purl.obolibrary.org/obo/go.owl")#, **kwargs)
            else:
                print("Using package release ontology")
                setOntology(owlFile=f"{os.path.dirname(os.path.abspath(__file__))}/data/go.owl")
        except Exception as e:
            print(f"Could not create ontology")
            print(e)
        
        self.tree = createGoTree(          ns       = ns,
                                  proteinList       = proteinList, 
                                  uniprotCollection = backgroundUniColl)

        proteinUniverseList = [ k for k in backgroundUniColl.keys() ]
        self.tree_universe = createGoTree(                ns = ns, 
                                           proteinList       = proteinUniverseList, 
                                           uniprotCollection = backgroundUniColl)

    @property
    def dimensions(self):
        return (self.tree.dimensions, self.tree_universe.dimensions)
    
    def walk(self, type=None):
        """Iterate over the tree, by default the experimental one"""
        if type is None:
            for node in self.tree.walk():
                yield node

    def computeORA(self, proteinIDList, rootPathway=None, ns="biological process", verbose=False):
        """Recursively computes the Fisher exact test of all pathways 
            based on corrresponding intersections between the following two categories
                    | Pa  | non_PA |
             -----------------------
                SA  |     |        |
              nonSA |     |        |
            Where,
                PA is the set of protein annotated in a pathway, extracted from the whole proteome tree
                SA is the set of proteins with altered experimental quantities: the <proteinIDList> parameter

        """
        pathWayXP = self.tree.getByName(ns)
        if rootPathway:
            pathWayXP = self.tree.getByName(rootPathway)
            if pathWayXP is None:
                print(f"No GO term named{None} was found in current ontology (maybe not populated)")
                return None
        
        pathWayBKG = self.tree_universe.getByName(ns)

        oraFisher, oraCDF = stat_utils.computeORA_BKG( pathWayXP,      # Pathway Node with all experimental proteins
                                                       proteinIDList,  # Set of protein w/ altered xptl qtty: observation counts
                                                       pathWayBKG,     # Available sets of protein in the pathways: hidden success counts [FROM PROTEOME]
                                                       verbose=verbose)
        return ORA_report(self.tree, oraFisher, oraCDF)
# Define rich view of stat results for notebook ?
class ORA_report:
    def __init__(self, xpTree, scoreF, scoreC):
        
        self.scoreF =  [ ORA_score(_) for _ in sorted(scoreF, key=lambda x:x[0]) ]
        self.scoreC =  [ ORA_score(_) for _ in sorted(scoreC, key=lambda x:x[0]) ]

    def __getitem__(self, x):
        return self.scoreF[x]

class ORA_score:
    def __init__(self, scoreElem):
        self.score, self.cNode = scoreElem
        self.proteins = self.cNode.getMembers()
    def __repr__(self):
        return f"{self.cNode.name} {self.score} {self.proteins}"

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