UNIVERSAL_TREES = {}

UNIVERSAL_VECTORS = {}

def clear():
    global UNIVERSAL_TREES, UNIVERSAL_VECTORS
    print(f"Clearing local stores content")

    UNIVERSAL_VECTORS = {}
    UNIVERSAL_TREES   = {}

def delTreeByTaxids(*taxids):
    miss = []
    for taxid in taxids:
        if not taxid in UNIVERSAL_TREES:
            miss.append(taxid)
        else:
            UNIVERSAL_TREES.pop(taxid)
    if miss:
        raise KeyError(f"{miss} to delete element tree not found in local store")

def delVectorsByTaxid(*taxids):
    miss = []
    for taxid in taxids:
        if not taxid in UNIVERSAL_VECTORS:
            miss.append(taxid)
        else:
            UNIVERSAL_VECTORS.pop(taxid)
    
    print(f"Following vector element were deleted {set(taxids) - set(miss)}")
    if miss:
        raise KeyError(f"{miss} to delete element vector not found in local store")

def storeTreeByTaxid(tree, taxid):
    print(f"localTreeStoring -->{taxid}")

    global UNIVERSAL_TREES
    if taxid in UNIVERSAL_TREES:
        raise KeyError(f"{taxid} tree already exists in local store")
        
    UNIVERSAL_TREES[taxid] = tree

def getUniversalTree(fullTreeKey):
    if not fullTreeKey in UNIVERSAL_TREES:
        raise KeyError(f"No tree names {fullTreeKey} found in local store")
    
    return UNIVERSAL_TREES[fullTreeKey]

def storeVector(vector, fullVectorKey):
    print(f"localVectorStoring -->{fullVectorKey}")

    global UNIVERSAL_VECTORS
    if fullVectorKey in UNIVERSAL_VECTORS:
        raise KeyError(f"{fullVectorKey} vector already exists in local store")
        
    UNIVERSAL_VECTORS[taxid] = tree

def getUniversalVector(taxid):
    if not taxid in UNIVERSAL_VECTORS:
        raise KeyError(f"No taxid {taxid} found in vector local store")
    
    return UNIVERSAL_VECTORS[taxid]

def getTaxidKeys():
    return list(UNIVERSAL_TREES.keys())

def listTreeKey(*args,**kwargs):
    return list(UNIVERSAL_TREES.keys())

def listVectorKey(*args,**kwargs):
    return list(UNIVERSAL_VECTORS.keys())
