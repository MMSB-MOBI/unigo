from . import redis
from . import local

CACHE_PKG=local
T_CACHE_TYPE=["local", "redis"]

def setCacheType(_type, **kwargs):
    if _type not in T_CACHE_TYPE:
        raise ValueError(f"{_type} is not registred cache type {T_CACHE_TYPE}")
    global CACHE_PKG
    CACHE_PKG = redis if _type == "redis" else local

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
    print(f"delTreeByTaxids :: {taxids}")
    return CACHE_PKG.delTreeByTaxids(taxids)

def getTaxidKeys():
    return CACHE_PKG.getTaxidKeys()

def getUniversalTree(taxid, raw=False):
    # if not local deserialize or deeper

    return CACHE_PKG.getUniversalTree(taxid, raw=raw)

def getUniversalVector(taxid):
    try:
        vec = CACHE_PKG.getUniversalVector(taxid)
        print(f"Found in vector cache {taxid}")
    except KeyError:
        try:
            tree = CACHE_PKG.getUniversalTree(taxid)
        except KeyError:
            raise KeyError(f"Vector error, No taxid {taxid} in stores")

        print(f"Creating vectors for {taxid}")
        vec = tree.vectorize()
        CACHE_PKG.storeVectorByTaxid(vec, taxid)
    
    return vec

def wipe():
    CACHE_PKG.wipe()

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