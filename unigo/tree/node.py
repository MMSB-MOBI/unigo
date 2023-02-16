import re, json, copy
from .heap import CoreHeap as kNodes 

NSGO = ["GO:0005575", "GO:0003674", "GO:0008150"]

class Node():

    def __init__(self, ID, name, oNode=None):
        self.ID = ID 
        self.name = name
        self.eTag =  [] # List of elements actually carrying the annotation ("tagged by the nodeName/annotation")
        self.leafCount =  0
        self.children =  []
        self.features = {}
        self.oNode = oNode
        self.isDAGelem = False
        self.is_a = [] # Used to deserialize from api
        self.background_frequency = None
        self.background_members =  [] # List of "members" elements of this node in the background proteome
        #self.heap = None

    def __deepcopy__(self, memo):
        # Deepcopy only the id attribute, then construct the new instance and map
        # the id() of the existing copy to the new instance in the memo dictionary
        memo[id(self)] = newself = self.__class__(copy.deepcopy(self.ID, memo), copy.deepcopy(self.name, memo), copy.deepcopy(self.oNode, memo))
        # Now that memo is populated with a hashable instance, copy the other attributes:
        newself.eTag = copy.deepcopy(self.eTag, memo)
        newself.background_members = copy.deepcopy(self.background_members, memo)
        newself.background_frequency = copy.deepcopy(self.background_frequency, memo)
        # Safe to deepcopy now, because backreferences to self will
        # be remapped to newself automatically
        newself.children = copy.deepcopy(self.children, memo)
        newself.features = copy.deepcopy(self.features, memo)
        newself.isDAGelem = copy.deepcopy(self.isDAGelem, memo)
        newself.oNode = copy.deepcopy(self.oNode, memo)
        #newself.heap = copy.deepcopy(self.heap, memo)

        return newself


    #def serial(self):


    def __hash__(self):
        return hash(self.ID)
    
    def __eq__(self, other):
        return hash(self) == hash(other)

# We may have to memo these 3   
    def as_DAG(self):
        return {
            'id' : self.ID,
            'children' : [ c.as_DAG() for c in self.children ]
        }

    def _as_DAG_strat(self, flat):
        if len(self.children) == 0 :
            return

        for c in self.children:
            if not c.ID in flat:
                flat[c.ID] = { "id" : c.ID, "parentIds" : set() }
            flat[c.ID]["parentIds"].add(self.ID)
        
        for c in self.children:
            c._as_DAG_strat(flat)
    
    def as_DAG_strat(self):
        #stHeap = kNodes()
        flat = { self.ID : { "id": self.ID, "parentIds" : [] } }
        #stHeap.add(self)

        self._as_DAG_strat(flat)
        
        return [ { "id" : v["id"], "parentIds" : list( v["parentIds"] )}  for k,v in flat.items() ]

    @property
    def pvalue(self):
        if 'Fisher' in self.features:
            return self.features['Fisher']
        return None
        

    def set(self, **kwargs):
        for k,v in kwargs.items():
            self.features[k] = v

    def __getattr__(self, key):
      #  print(key)
    #   if key not in self.features:
        if key == "ID":
            print(dir(self))
            return self.__getattribute__(key)

        if key == "features":
            raise AttributeError(key)
        try:
            return self.features[key]
        except KeyError:
            raise AttributeError(key)

        #return self.features[key]
    
    def __dir__(self):
        return super().__dir__() + [str(k) for k in self.features.keys()]
    
    def hasChild(self, ID):
    # print(node)
        for child in self.children:
            if child.ID == ID:
                return child
        return None
    
   
    def traverse(self):
        yield self
        for c in self.children:
            yield from c.traverse()
     
    # Memoized version of traverse
    def walk(self, **kwargs):
        wHeap = kNodes()
        return self._walk(wHeap, **kwargs)

    def _walk(self, wHeap, **kwargs):
        if self.isDAGelem and self in wHeap:            
            return        
        wHeap.add(self) 
        if 'mustContain' in kwargs:
            if not set(kwargs["mustContain"]) & set(self.getMembers()) :
                #print(f"-->Early exit at empty node {self.name}")
                return

        yield self
        for c in self.children:
            yield from c._walk(wHeap, **kwargs)

    def _as_newick(self):
        #_self = node['name'].replace(" ", "_")
        _self = self.name.replace(',',' ').replace('(', '[').replace(')', ']').replace('\'', '_').replace(':', '_')
        
        if len(self.children) == 0:
            return _self

        return '(' + ','.join([ c._as_newick() for c in self.children ]) + ')' + _self
            
    def getByName(self, name):       
        regExp = name.replace(' ', '.').replace('[', '.').replace(']', '.').replace('_', '.').replace(')', '\)').replace('(', '\(')
        regExp = regExp.replace('+', '.')
        if re.search("^" + regExp + "$", self.name):
            return self
        for c in self.children:
            v = c.getByName(name)
            if v:
                return v
        return None

    def getByID(self, ID):       
        if self.ID == ID:
            return self
        for c in self.children:
            v = c.getByID(ID)
            if v:
                return v
        return None


# Memoization implies the leaves of the subtree are already registred
# API to perform search at current node, we clear the heap

    def getMembers(self, nr=False):
        getMemberHeap = kNodes()
        getMemberHeap.add(self)

        buff = copy.copy(self.eTag)
        for n in self.children:
            buff += n._getMembers(getMemberHeap)
        
        return buff if not nr else list(set(buff))

    def _getMembers(self, _heap):
        if self.isDAGelem and self in _heap:
            return []        
        _heap.add(self)
        buff = copy.copy(self.eTag)
        for n in self.children:
            buff += n._getMembers(_heap)
        
        return buff


    #def _mayCollapseNode(self, fnPredicate, _ctHeap):
    #""" Apply a predicate to collapse a node
    #    Collapsing a node makes 
    #"""


    def _collapseNode(self, _ctHeap):
        #if self.isDAGelem and self in _ctHeap:
        #    print(f"{self.name} already visited") ## This is buggy
        #    return self
        
        #print(f"Looking in {self.name}")
        #for n in self.children:
        #    print(f"{self.name} -> {n.name}")
        _ctHeap.add(self)
        #print(f"collapsing {self.name}")
        if len(self.children) == 0: # Leave or a node carrying actual protein, return it
            return self
           
        if len(self.children) == 1 and len(self.eTag) == 0 and self.ID not in NSGO:
            #print(f"Skipping {self.name}")
            return self.children[0]._collapseNode(_ctHeap)

        #print(f"Keeping {self.name} children={len(self.children)} eTag{len(self.eTag)}")
        self.children = [ n._collapseNode(_ctHeap) for n in self.children ]
        #print(self.children)
        
        return self

    def __repr__(self):
        d = { k : v for k,v in self.__dict__.items()  if not k == "children" }
        d["children"] = [ c.name for c in self.children ]
        return str(d)

    def __str__(self):
        return self.__repr__()

    def _mayDrop(self, predicate, _noDropHeap):
        #print(f"Testing {self.name}")
        if not predicate(self):
        #   print(f"{self.name} is droped !")
            return None

        if self in _noDropHeap:
            return self
        _noDropHeap.add(self)

        self.children = [ c for c in self.children if c.mayDrop(predicate, _noDropHeap) ]

        return self

    def mayDrop(self, predicate, pNode, _noDropHeap):
        # Self failed add nothing to stack
        if not predicate(self):        
            return
        # add self to stack
        _noDropHeap.add(self, pNode)

        # pass on to children
        for c in self.children:
            c.mayDrop(predicate, pNode, _noDropHeap)

        return

    def _leafCountUpdate(self, _lcHeap):
      #  print(f"lcu {self.name} {len(_lcHeap)}")
        if self.isDAGelem and self in _lcHeap:
            #print(f"{self.name} already updated")
            return self
        _lcHeap.add(self)

        self.leafCount = len(self.getMembers())
        #print(f"{self.name} leafCount is {self.leafCount}")

        for c in self.children:
            c._leafCountUpdate(_lcHeap)
