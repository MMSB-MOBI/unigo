from pyproteinsext import ontology
import uuid, os, pickle
import os.path
from .node import Node as createNode
from . import heap 
from owlready2 import get_ontology

GO_ONTOLOGY = None
GO_DICT = {}

enumNS = {
            'biological process' : 'GO:0008150',
            'molecular function' : 'GO:0003674',
            'cellular component' : 'GO:0005575'
        }

enumNSkeys = { # Usefull for redis key indexation
            'biological process' : 'P',
            'molecular function' : 'F',
            'cellular component' : 'C'
        }

enumNSkeysRevert = { # Usefull for http results wraping
            v:k for k, v in enumNSkeys.items()
        }


def assertAndCoherceValidNamespace(k):
    if not k in enumNSkeys:
        raise KeyError(f"{k} is not a valid GO namespace")
    return enumNSkeys[k]

def setOntologyDict(owlFile):
    global GO_DICT
    GO_DICT = {}
    ontology = get_ontology(owlFile).load()
    for node in ontology.classes():
        for node_id in set(node.id):
            if node_id not in GO_DICT: 
                GO_DICT[node_id] = node
            else:
                print("WARN", node_id, "already stored")
                
        if node.hasAlternativeId:
            for node_id in set(node.hasAlternativeId):
                if node_id not in GO_DICT:
                    GO_DICT[node_id] = node
                else:
                    print("WARN", node_id, "already stored")

        
def setOntology(owlFile=None, url=None):
    global GO_ONTOLOGY
    if not owlFile and not url:
        raise ValueError("file or url named parameters required")
    if owlFile:
        GO_ONTOLOGY = OntologyGO(owlFile=owlFile)
    else:
        GO_ONTOLOGY = OntologyGO(url=url)
    if GO_ONTOLOGY is None:
        raise ValueError(f"Could not set GO ontology from {owlFile} or {url}")

def ontologyDump(name, location):
    name = str( uuid.uuid4() ) if not name else name
    location = os.getcwd() if not location else location
    fName = location + '/' + name + '.sqlite3'
    print (f"Trying to open {fName}")
    ontology.default_world.set_backend(filename = fName)
    ontology.default_world.save()



# DAG implemnetaion requires memoizing visited node
# -> Any external accessor must be at the tree level, not at the node level avoid memoizing conflict
# Memoizing can be with a node heap or a node attribute
# We try the second one

# Check and raise Error at the 1st encountered cycle in
# the subtree of provided node 
import copy
def _checkCycle(rootNode, cNode, _path):
    path = copy.deepcopy(_path)
    path.append(cNode.ID)
    for child in cNode.children:
        if child.ID == rootNode.ID:
            pathString = '->'.join(path) + '->' + child.ID
            print (pathString)
            raise ValueError(f"{rootNode.name} has a cycle {pathString}")

def checkCycle(rootNode, verbose=False):
    if verbose:
        print(f"Checking {rootNode.ID}")

    for child in rootNode.children:
        _checkCycle(rootNode, child, [rootNode.ID])
    
    for child in rootNode.children:
        checkCycle(child)


class OntologyGO():
    def __init__(self, owlFile=None, url=None):
        if not owlFile is None:
            self.onto = ontology.Ontology(owlFile=owlFile)
        elif not url is None:
            self.onto = ontology.Ontology(url=url)
        else :
            raise("required parameter owlFile or url")
        
    def getLineage(self, goID):
        lin = self.onto._getLineage(goID)
        if not lin:
            return None
        return [ (t.id[0], t.label[0]) for t in lin[:-1] ]


def ascend(cNode, nodeSet, rootSet, strict):#;tree):
    """ Recursively look for the 1st ascendant node w/out parent
        
    """

    #print(f"-->{cNode.ID}")
    for _p in cNode.oNode.is_a:
        # non ascendant is_a element
        #http://www.ontobee.org/ontology/RO?iri=http://purl.obolibrary.org/obo/BFO_0000051
        if strict and not str(_p).startswith('obo.GO') and not str(_p) == 'owl.Thing':
            continue

        if str(_p).startswith('obo.BFO_0000051'):
            continue
        
        try : # trying to dereference Restriction wrapper
            p = _p.value
        except AttributeError:
            p = _p

        # cNode is terminal node
        if str(p) == "owl.Thing" :
            if len(cNode.oNode.is_a)> 1:
                raise TypeError(f"Terminal node having multiple parent: {cNode}")
            #print(f"Adding {cNode} as root")
            rootSet.add(cNode)
            return
        #print(f"{p}")
        # the current pNode is already registred, we eventually register cNode as one of its children
        # But won't keep on ascending
        pNode =  createNode(p.id[0], p.label[0], p)
        toStop = True  if pNode in nodeSet else False

        # register current parent node
        pNode = nodeSet.add( pNode )
        pNode.isDAGelem = True
        #pNode.heap = tree.nodeHeap
        
        if not pNode.hasChild(cNode.ID):
            pNode.children.append(cNode)
        
        if not toStop:
            ascend(pNode, nodeSet, rootSet, strict)#, tree)

def collapseTree(_root):
    root = copy.deepcopy(_root)
    ctHeap = heap.CoreHeap()
    root.children = [ n._collapseNode(ctHeap) for n in root.children ]
    return root
    
def insertLineage(root, lineage, eName):
    cNode = root
    for (goID, name) in reversed(lineage):
        mNode = cNode.hasChild(goID)
        if not mNode:
            #print(f"{goID},{name} NOT found under {cNode['ID']},{cNode['name']}")
            mNode = createNode(goID, name, None)
            cNode.children.append(mNode)
        else :
            #print(f"{goID},{name} found under {cNode['ID']},{cNode['name']}")
            pass
        cNode = mNode
        cNode.leafCount += 1
    cNode.eTag.append(eName)
    cNode.leafCount -= 1

def deserializeGoTree(fPickle, owlFile):
    global GO_ONTOLOGY
    if not GO_ONTOLOGY:
        setOntology(owlFile)
    
    fp = open(fPickle, 'rb')
    _self = pickle.load(fp)
    fp.close()
    for n in _self.walk():
        if n.oNode is None:
            print(f"Root found skipping::{n}\n")
            continue
        if isinstance(n.oNode, str):
            termID = n.oNode#.replace('obo.', '').replace('_', ':')
            termObj = GO_ONTOLOGY.onto.onto.search_one(id=termID)

            if termObj is None:
                raise TypeError(f"No GO Term matching::{termID}::From::{n}\n")
            n.oNode =  termObj

    return _self

def createGoTree(ns = None, protein_iterator = None, collapse = True, strict = True):
    if  ns is None:
        raise ValueError("Specify a namespace \"ns\"")
    if protein_iterator is None:
        raise ValueError("Specify a protein iterator \"protein_iterator\"")

    xpGoTree = AnnotationTree(ns, collapse)
    print(f"Blueprint xpGoTree {ns} extracted")
    xpGoTree.extract(protein_iterator, strict)
    print(f"xpGoTree {ns} filtered for supplied uniprot entries, indexing by uniprot_ids")
    xpGoTree._index()
    return xpGoTree

def load(baseData):
    """Deserialize tree as dict structure fetched from api
        parameter is a shallow dict structure wich is assumed to 
        feature a single value with an empty "is_a" array ie a single root element
    {
        "0000":{"children":["GO:0008150"],"eTag":[],"is_a":[],"name":"root"},
        "GO:0008150": {"children":["GO:0051179","GO:0009987","GO:0008152","GO:0050896"],"eTag":[],"is_a":["0000"],"name":"biological_process"},
        ...
    }
    """
    
    flatNodeData = spawn(baseData)
    root = wire(flatNodeData, baseData)
    maybeNS = root.children[0].name.replace('_', ' ')
    
    t = AnnotationTree(maybeNS, collapse=True)
    t.root = root
    
    t.isDAG = True
    t.collapsable = True # Arbitrary
    t.leafCountUpdate()
    #TO DO ?
    # print("Indexing uniprot identifiers")
    #t._index()
    return t
        
def spawn(strData):
    """Utility deserialize function
        Generate a shallow dictionary of Node objects
    """
    nodeData = {}
    for nKey, nDatum in strData.items():
        nodeData[nKey] = createNode(nKey, nDatum['name'])
        nodeData[nKey].eTag = nDatum['eTag']
        nodeData[nKey].isDAGelem = nDatum['isDAGelem']
        nodeData[nKey].background_frequency = nDatum['background_frequency']
        nodeData[nKey].background_members = nDatum['background_members']
    #print(f"spawn Nodedata has {len(list(nodeData.keys()))} elem")
    
    return nodeData

def wire(nodeData, strData, mayDropOccur=False):
    """Utility deserialize function
        fill Node reference across provided dictionary on Node objects
        using Node identifiers stored in 2nd paramter
        return the root of corresponding tree
    """
    root=None
    for nKey, nObj in nodeData.items():
        nRefStr = strData[nKey]
        if len(nRefStr['is_a']) == 0:
            if not root is None:
                raise KeyError("Multiple roots")
            root = nObj
        for childID in nRefStr['children']: # Children missing reference are allowed in the context of drop based tree reconstruction
            if not childID in nodeData:
                if mayDropOccur:
                    continue
                else:
                    raise KeyError(f"while editing children reference, {childID} not found in shallow node dictionary")
            nObj.children.append(nodeData[childID])

        for parentID in nRefStr['is_a']: # parent missing reference are forbidden for now
            nObj.is_a.append(nodeData[parentID])
    
    return root

   



class AnnotationTree():
    def __init__(self, annotType,  collapse=False):
        
        if annotType not in enumNS:
            raise KeyError (f"annotation type \"{annotType}\" is not allowed ({enumNS}) {{ {enumNS.keys()} }}")

        self.collapsable = collapse
        self.isDAG = False
        #self.nodeHeap = heap.CoreHeap()
        self.root = createNode('0000', 'root') 
        #self.root.heap = self.nodeHeap
        self.NS = (annotType, enumNS[annotType])
        self.index_by_uniprotid = None
    
    def compute_background_frequency(self):
        nb_all_proteins = self.dimensions[3]
        for n in self.walk():
            nb_node_proteins = len(set(n.getMembers()))
            frq = nb_node_proteins / nb_all_proteins
            n.background_frequency = frq


    def vectorize(self):
        """Create a dictionary representation of passed tree content, getting rid of tree topology 
        in exchanged for a vectorizable data structure of node terms lited under the "elements" list.
        All dictionary key,value are strings or integer for easy serialization.
        The returned dictionary layout is the following:

        { 
            "registry" : [uniprotID,],
            "terms" : { "GOid" : [ elements.index, ... ], }
        } 

        Parameters
        ----------
        aTree: a tree.AnnotationTree object.

        Returns
        -------
        A dictionary of tree nodes
        """
        vHeap = heap.VHeap()
        for n in self.walk():
            vHeap.add(n, None)
        return vHeap.asDict

    # Serializing without any ONTOLOGY references
    def f_serialize(self):
        fHeap = heap.FHeap()
        cnt = 0
        for cNode in self.walk():
            fHeap.add(cNode, None)
            for node in cNode.children:
                fHeap.add(node, cNode)
            cnt += 1
        print(f"\n#################\nAnnotation tree serialization")
        print(f"NS:{self.NS}")
        print(f"{cnt} nodes traversal")
        print(f"heap dimensions:: nodes:{fHeap.dimensions[0]} childrenRefCount :{fHeap.dimensions[1]}")
        print(f"#################\n")   
        return fHeap
    
    # Serialization requires pickling of object instance
    # Ontology related term will be saved as string
    # The ontology object will have to be reimported for deserialization
    def makePickable(self):
        _self = copy.deepcopy(self)
        for n in _self.walk():
            n.oNode = str(n.oNode).replace('obo.', '').replace('_', ':') if not n.oNode is None else None
        return _self

    def extract(self, protein_iterator, strict):
        self.read_DAG(protein_iterator, strict)
        
    def _index(self):
        if self.index_by_uniprotid is None:
            self.index_by_uniprotid = {}
        
        for cNode in self.walk():
            for uniprotid in cNode.getMembers():
                if not uniprotid in self.index_by_uniprotid:
                    self.index_by_uniprotid[uniprotid] = []
                self.index_by_uniprotid[uniprotid].append(cNode)
    

        #walk(mustContains=)
        
    def read_DAG(self, uniprot_iterator, strict):  
        print('read dag')    
        """ Cross GO Ontology with supplied uniprot_iterator
            to create the minimal GO DAG containg all GO terms featured by uniprot collection
        """
        self.isDAG = True
        global GO_ONTOLOGY
        if GO_ONTOLOGY is None:
            print("Please set GO_ONTOLOGY")
            return
      
        ontologyNode = GO_ONTOLOGY.onto.onto.search_one(id=self.NS[1])
        if not ontologyNode:
            raise ValueError(f"id {enumNS[self.NS[1]]} not found in ontology")
      #  self.root.children.append( Node(enumNS[annotType], annotType, oNode=ontologyNode) )

        nodeSet = heap.CoreHeap()
        rootSet = heap.CoreHeap()
        def setSentinelChar():
            """ Returns the letter used by GO to prefix its term depending on namespace"""
            if self.NS[0] == 'biological process':
                return 'P'
            elif self.NS[0] == 'molecular function':
                return 'F'
            return 'C'
        goNSasChar = setSentinelChar()
        disc = 0
        i = 0
        for prot in uniprot_iterator:
            goTerms = prot.go
            uniID = prot.id
            bp = []
            for goTerm in goTerms:
                if goTerm.term.startswith(f"{goNSasChar}:"):
                    bp.append(goTerm.id)
            if not bp:
                disc += 1
                #print(f"Added {p} provided not GO annnotation (current NS is {self.NS[0]})")
            for term in bp:
                cLeaf = GO_ONTOLOGY.onto.onto.search_one(id=term)
                if not cLeaf:
                    cLeaf = GO_ONTOLOGY.onto.onto.search_one(hasAlternativeId=term)
                    if not cLeaf:                       
                        print("Warning: " + term + " not found in "+\
                             self.NS[0] + ", plz check for its deprecation "+\
                             "at " + "https://www.ebi.ac.uk/QuickGO/term/" + term)
                        continue
                #print(f"adding {term}")
                #print(f"with// createNode({cLeaf.id[0]}, {cLeaf.label[0]}, {cLeaf}")
                # Add a new node to set of fetch existing one
                bottomNode = nodeSet.add( createNode(cLeaf.id[0], cLeaf.label[0], cLeaf) )
                bottomNode.eTag.append(uniID)
                bottomNode.isDAGelem = True
                #bottomNode.heap = self.nodeHeap
                #print(f"rolling up for {bottomNode.ID}")
                #print(f"rolling up {term}")
                ascend(bottomNode, nodeSet, rootSet, strict)#, self)
                #print(f"S1a stop {term}")
            #print(f"{uniID} done")
        #if len(rootSet) > 1:
        #    raise ValueError(f"Too many roots ({len(rootSet)}) {list(rootSet)}")

        

        for n in rootSet:
            if n.ID == self.NS[1]:
                self.root.children.append(n)

        if self.collapsable:
            print("Applying true path collapsing")
            self.root = collapseTree(self.root)
            #self.nodeHeap = self.root.heap
            
        n, ln, l, p = self.dimensions
        print(f"{n} GO terms, {ln} children_links, {l} leaves, {p} proteins ({disc} discarded)")     

    def read(self,uniprotIDList, uniprotCollection):      
        global GO_ONTOLOGY
        if GO_ONTOLOGY is None:
            print("Please set GO_ONTOLOGY")
            return
        
        self.root.children.append( createNode(self.NS[1], self.NS[0]) )
               
        i=0
        for p in uniprotIDList:
            u = uniprotCollection.get(p)
            bp = list((u.GO[self.NS[0]]).keys())
            if bp:
                i += 1          
            for term in bp:
                lineage = [(term, u.GO[self.NS[0]][term])] + GO_ONTOLOGY.getLineage(term)
                insertLineage(self.root, lineage, p)
                
        print(f"Annotation {self.NS[0]} extracted from {i} / {len(uniprotCollection)} uniprot entries parsed")
        
        if self.collapsable:
            print("Applying true path collapsing")
            self.root = collapseTree(self.root)
          #  self.nodeHeap = self.root.heap
    
    @property 
    def dimensions(self):
        nNodes = 0
        nLinks = 0
        for n in self.root.walk():
            nNodes += 1
            nLinks += len(n.children)
        leafTotal = self.root.getMembers()
        leafTotal_nr = set(leafTotal)

        return nNodes, nLinks, len(leafTotal), len(leafTotal_nr)

    @property
    def proteins(self):
        return set(self.root.getMembers())

    def walk_no_root(self):
        skip_root = False
        for node in self.root.walk():
            if not skip_root:
                skip_root = True
                continue
            yield node
    
    def walk(self, **kwargs):       
        return self.root.walk(**kwargs) 
    
    def as_json(self): # Should be collapsible
        return json.dumps(self.root)

    def as_newick(self, collapse=True):
        return '(' + self.root._as_newick() + ');'

    def getByName(self, name):
        n = self.root.getByName(name)
        if not n:
            raise KeyError(f"No node named \"{name}\" in current tree")
        return n

    def _getByID(self, ID):

        n = self.root.getByID(ID)
        if not n:
            raise KeyError(f"No node w/ ID \"{ID}\" in current tree")
        return n

    def collapse(self):
        t = AnnotationTree(self.NS[0])
        t.isDAG = self.isDAG
      

        t.root = collapseTree(self.root)
       # t.nodeHeap = t.root.heap
        t.collapsable = True
        return t

    def getMembersByName(self, name, nr=True):
        node = self.root.getByName(name)
        m = node.getMembers()
        if nr:
            s = set(m)
            return list(set(m)) 

        return list(m)
       
    def getMembersByID(self, name, nr=True):
        node = self.root.getByID(name)
        #print(node)
        m = node.getMembers()
        if nr:
            return list(set(m))
        
        return list(m)

    def getMembers(self):
        return self.root.getMembers(nr=True)

    def getDetailedMembersFromParentID(self, parent_id):
        def addToMembers(node):
            if node.ID not in browsed:
                for member in node.eTag:
                    if member not in members:
                        members[member] = []
                    members[member].append((node.ID, node.name))
                browsed.append(node.ID)
            if node.children:
                for child in node.children:
                    addToMembers(child)

        browsed = []
        members = {}
        parent_node = self.getByID(parent_id)
        addToMembers(parent_node)
        return members

    def newRoot(self,  name=None, ID=None):
        if not name and not ID:
            raise ValueError("Provide name or ID")

        t = AnnotationTree(self.NS[0])
        t.isDAG = self.isDAG
        t.collapsable = self.collapsable

        field = "name" if not name is None else "ID"
        value = name if not name is None else ID
        fn = self.getByName if not name is None else self.getByID

        nr = fn(value)

        if not nr:
            raise KeyError(f"No node with {field} {value} in mother tree")
        t.root.children.append(nr)
        
        return t

# A node may effectively be tested more than once
# Given the operation cost associated with False predicate (terminate&return None), heap is unecessary
# Given the operation cost associated with True predicate we effectively test the below tree once again
    def _drop(self, predicate, noCollapse=False):
        t = AnnotationTree(self.NS[0])
        t.isDAG = self.isDAG
        t.collapsable = self.collapsable
 
        noDropHeap = heap.CoreHeap() # To store success and avoid restest subtree

        t.root = copy.deepcopy(self.root)
        t.root.children = [ c for c in t.root.children if c._mayDrop(predicate, noDropHeap) ]

        if not noCollapse:
            t = t.collapse()
        t.leafCountUpdate()

        return t
    
    def drop(self, predicate, noCollapse=False, noLeafCountUpdate=False):
        """Returns a subtree with all nodes matching parameter predicate function
        """
        noDropHeap = heap.FHeap() # To store success and avoid restest subtree
        noDropHeap.add(self.root, None)
        # start recursive
        for _ in self.root.children:
            _.mayDrop(predicate, self.root, noDropHeap)
        
        # Build results tree from the flat dictionary of surviving nodes
        baseData = noDropHeap.asDict
        flatNodeData = spawn(baseData)
        root = wire(flatNodeData, baseData, mayDropOccur=True)
        
        t = AnnotationTree(self.NS[0])
        t.isDAG = self.isDAG
        t.collapsable = self.collapsable
        t.root = root
        if not noCollapse:
            t = t.collapse()
        if not noLeafCountUpdate:
            t.leafCountUpdate()

        return t

    def leafCountUpdate(self):
        lcHeap = heap.CoreHeap()

     #   if self.isDAG:
     #       print(f"clearing heap")
     #       self.nodeHeap.clear()
     #       print(f"heap size: {len(self.nodeHeap)}")
        for c in self.root.children:
            c._leafCountUpdate(lcHeap)

    def as_DAG_strat(self):
        return self.root.as_DAG_strat()
