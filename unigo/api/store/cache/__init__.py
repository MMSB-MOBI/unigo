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
    #print(f"Set cache to {_type}")

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

def delVectorsByTaxid(taxids, delAllRelated=True):
    CACHE_PKG.delVectorsByTaxid(taxids)
    if CACHE_SYMBOL == 'redis' and delAllRelated:
        CACHE_PKG.delCulledByTaxids(taxids)

def getTaxidKeys():
    return CACHE_PKG.getTaxidKeys()

def getUniversalTrees(taxid, raw=False):
    # if not local deserialize or deeper

    return CACHE_PKG.getUniversalTree3NS(taxid, raw=raw)

def getUniversalVectors(taxid):
    try:
        vectors = CACHE_PKG.getUniversalVector3NS(taxid)
        print(f"{taxid} Universal vector in cache")
    except KeyError:
        try:
            tree3NS = CACHE_PKG.getUniversalTree3NS(taxid)
        except KeyError:
            raise KeyError(f"Vector error, No taxid {taxid} in stores")

        print(f"Building 3NS {taxid} Universal vector")
        vectors = {}
        for tree in tree3NS:
            _ = tree.vectorize()
            vectors[ tree.ns ] = _
            CACHE_PKG.storeVectorByTaxid(_, taxid)
    
    return vectors

### FOLLOWING TWO SHOULD BE REMOVED ###
def _getUniversalTree(taxid, raw=False):
    # if not local deserialize or deeper

    return CACHE_PKG.getUniversalTree(taxid, raw=raw)

def _getUniversalVector(taxid):
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

def deleteTaxids(taxids):
    return CACHE_PKG.deleteTaxids(taxids)

def status():
    if CACHE_SYMBOL == 'redis':
        return len(listTrees()), len(listVectors()), len(listCulled())

def listTrees():
    if CACHE_SYMBOL == 'redis':
        return [ _ for _ in CACHE_PKG.listTreeKey(prefix=False) ]

def listVectors():
    if CACHE_SYMBOL == 'redis':
        return [ _ for _ in CACHE_PKG.listVectorKey(prefix=False) ]
    else:
        raise TypeError("YOU SHOULD IMPLEMENT LOCAL KEYS ITER")

def listCulled():
    if CACHE_SYMBOL == 'redis':
        return [ _ for _ in CACHE_PKG.listCulledVectorKey(prefix=False) ]
    else:
        raise TypeError("YOU SHOULD IMPLEMENT LOCAL KEYS ITER")


def getCulledVectors(taxid, cmin, cmax, fmax):
    if CACHE_SYMBOL == 'redis':
        _ = CACHE_PKG.getCulledVector3NS(taxid, cmin, cmax, fmax)      
    else:
        raise TypeError("YOU SHOULD IMPLEMENT LOCAL getCulledVector")
    return _

def storeCulledVector(vector, taxid, ns, cmin, cmax, fmax):
    if CACHE_SYMBOL == 'redis':
        _ = CACHE_PKG.storeCulledVector(vector, taxid, ns, cmin, cmax, fmax)
    else:
        raise TypeError("YOU SHOULD IMPLEMENT LOCAL storeCulledVector")
    return _   

def listMissUniversalVector():
    _treeID   = set( listTrees() )
    _vectorID = set( listVectors() )
    #print(f"{_treeID}  -----  {_vectorID}")
    return list(_treeID - _vectorID)

def buildUniversalVector():
    taxid_ns_treeKeys = listMissUniversalVector()
    for taxid_ns_treeKey in taxid_ns_treeKeys:
        tree = CACHE_PKG.getUniversalTree(taxid_ns_treeKey)
        print(f"Build vector for {taxid_ns_treeKey} from {tree}")
        CACHE_PKG.storeVector(tree, key=taxid_ns_treeKey)
