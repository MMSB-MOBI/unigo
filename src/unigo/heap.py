class CoreHeap():
    def __init__(self):
        self.data = {}

    def __len__(self):
        return len(list(self.data.keys()))

    def  __contains__(self, node):
        return node in self.data

    def add(self, node):
        if not node in self:
            self.data[node] = node
        return self.data[node]

    def clear(self):
        self.data = {}

    @property
    def asDict(self):
        return self.data

    @property
    def dimensions(self):
        lnkNum = 0
        for k,v in self.data.items():
            lnkNum += len(v['children']) 
        return len(self), lnkNum
    
    def __iter__(self):
        for k, v in self.data.items():
            yield v

    def items(self):
        for k, v in self.data.items():
            yield (k,v)       

class FHeap(CoreHeap):
    def __init__(self):
        super().__init__()
    def add(self, cNode, pNode):
        if not cNode.ID in self.data:
            self.data[cNode.ID] = {
                'name' : cNode.name,
                'eTag' : cNode.eTag,
                'children' : [ _.ID for _ in cNode.children ],
                'is_a' : [],
                'isDAGelem' : cNode.isDAGelem
            }
        if not pNode is None:
            self.data[cNode.ID]['is_a'].append(pNode.ID) # Do we check unicity?

class VHeap(CoreHeap):
    def __init__(self):
        super().__init__()
        self.totalElements = set()
    def add(self, cNode, pNode):
        if not cNode.ID in self.data:
            _  = cNode.getMembers(nr=True)
            self.data[cNode.ID] = {
                'name' : cNode.name,
                'elements' : _,
            }
            self.totalElements = self.totalElements | set(_)
        if not pNode is None:
            self.data[cNode.ID]['is_a'].append(pNode.ID) # Do we check unicity?
    
    @property
    def asDict(self):
        registry = list(self.totalElements)

        _ = { "registry"  : registry,
              "terms"     : {}
            }
        
        for goID, goTerm in self.items():
            _["terms"][goID] = {
                "name" : goTerm['name'],
                "elements" : [ registry.index(e) for e\
                                in  goTerm["elements"] ],
                "freq"   : float( len(goTerm["elements"])/ len(registry) )  
            }
        return _
