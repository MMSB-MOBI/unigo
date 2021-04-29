from . import redis
from . import local

CACHE_PKG=local
T_CACHE_TYPE=["local", "redis"]
CACHE_SYMBOL="local"

def setCacheType(_type, **kwargs):
    if _type not in T_CACHE_TYPE:
        raise ValueError(f"{_type} is not registred cache type {T_CACHE_TYPE}")
    global CACHE_PKG, CACHE_SYMBOL
    CACHE_PKG = redis if _type == "redis" else local
    CACHE_SYMBOL = _type
    print(f"Set cache to {_type}")

    if _type == 'redis':
        p = {}
        if 'rp' in kwargs:
            p['port'] = kwargs['rp']
        if 'rh' in kwargs:
            p['host'] = kwargs['rh']
        
        CACHE_PKG.setDatabaseParameters(**p)

def storeTreeByTaxid(tree, taxid):
    return CACHE_PKG.storeTreeByTaxid(tree, taxid)

def delTreeByTaxids(taxids):
    return CACHE_PKG.delTreeByTaxids(taxids)

def getTaxidKeys():
    return CACHE_PKG.getTaxidKeys()

def getUniversalTree(taxid, raw=False):
    # if not local deserialize or deeper

    return CACHE_PKG.getUniversalTree(taxid, raw=raw)

def getUniversalVector(taxid):
    try:
        vec = CACHE_PKG.getUniversalVector(taxid)
        print(f"{taxid} Universal vector in cache")
    except KeyError:
        try:
            tree = CACHE_PKG.getUniversalTree(taxid)
        except KeyError:
            raise KeyError(f"Vector error, No taxid {taxid} in stores")

        print(f"Building {taxid} Universal vector")
        vec = tree.vectorize()
        CACHE_PKG.storeVectorByTaxid(vec, taxid)
    
    return vec

def clear():
    CACHE_PKG.clear()

def status():
    if CACHE_SYMBOL == 'redis':
        return len(listTrees()), len(listVectors())

def listTrees():
    if CACHE_SYMBOL == 'redis':
        return [ _ for _ in CACHE_PKG.listTreeKey(prefix=False) ]

def listVectors():
    if CACHE_SYMBOL == 'redis':
        return [ _ for _ in CACHE_PKG.listVectorKey(prefix=False) ]
    else:
        raise TypeError("YOU SHOULD IMPLEMENT LOCAL KEYS ITER")

def getCulledVector(taxid, cmin, cmax, fmax):
    if CACHE_SYMBOL == 'redis':
        _ = CACHE_PKG.getCulledVector(taxid, cmin, cmax, fmax)
    else:
        raise TypeError("YOU SHOULD IMPLEMENT LOCAL getCulledVector")
    return _

def storeCulledVector(vector, taxid, cmin, cmax, fmax):
    if CACHE_SYMBOL == 'redis':
        _ = CACHE_PKG.storeCulledVector(vector, taxid, cmin, cmax, fmax)
    else:
        raise TypeError("YOU SHOULD IMPLEMENT LOCAL storeCulledVector")
    return _   

def listMissUniversalVector():
    _treeID   = set( listTrees() )
    _vectorID = set( listVectors() )
    return list(_treeID - _vectorID)

def buildUniversalVector():
    print("Running unBuildtUniversalVectorIter")
    _ = listMissUniversalVector()
    print(f"Set of trees to vectorize {_}")
    for bKey in _:
        print(f"-->{type(bKey)}")
        print(f"Build vector for {bKey}")
        tree = CACHE_PKG.getUniversalTree(bKey)
        CACHE_PKG.storeVectorByTaxid(tree, bKey)

"""
def storeTreeByTaxid(tree, taxid):
    return FN_HANDLER['storeTreeByTaxid'](tree, taxid)


def taxidKeys():
    return FN_HANDLER['getTaxidKeys']

def getUniversalTree(taxid):
    # if not local deserialize or deeper

    return FN_HANDLER['getUniversalTree'](taxid)

# Can be called for construction of vector banking
def getUniversalVector(taxid):
    try:
        vec = FN_HANDLER['getUniversalVector'](taxid)
        print(f"Found in vector cache {taxid}")
    except KeyError:
        try:
            tree = FN_HANDLER['getUniversalTree'](taxid)
        except KeyError:
            raise KeyError(f"Vector errro, No taxid {taxid} in stores")

        print(f"Creating vectors for {taxid}")
        vec = tree.vectorize()
        FN_HANDLER['storeVectorByTaxid'](vec, taxid)
    
    return vec
"""