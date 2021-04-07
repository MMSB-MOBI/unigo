UNIVERSAL_TREES = {}

UNIVERSAL_VECTORS = {}

def wipe():
    global UNIVERSAL_TREES, UNIVERSAL_VECTORS
    print(f"Wiping local stores content")

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
        raise KeyError(f"{miss} to delete elements not found in local store")

def storeTreeByTaxid(tree, taxid):
    print(f"localTreeStoring -->{taxid}")

    global UNIVERSAL_TREES
    if taxid in UNIVERSAL_TREES:
        raise KeyError(f"{taxid} tree already exists in local store")
        
    UNIVERSAL_TREES[taxid] = tree

def getUniversalTree(taxid):
    if not taxid in UNIVERSAL_TREES:
        raise KeyError(f"No taxid {taxid} found in local store")
    
    return UNIVERSAL_TREES[taxid]

def storeVectorByTaxid(vector, taxid):
    print(f"localVectorStoring -->{taxid}")

    global UNIVERSAL_VECTORS
    if taxid in UNIVERSAL_VECTORS:
        raise KeyError(f"{taxid} vector already exists in local store")
        
    UNIVERSAL_VECTORS[taxid] = tree

def getUniversalVector(taxid):
    if not taxid in UNIVERSAL_VECTORS:
        raise KeyError(f"No taxid {taxid} found in vector local store")
    
    return UNIVERSAL_VECTORS[taxid]

def getTaxidKeys():
    return list(UNIVERSAL_TREES.keys())
