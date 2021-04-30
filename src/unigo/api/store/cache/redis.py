from pyrediscore import connect, store, get, delete, listKey, setDatabaseParameters
from ...data_objects import loadUnivGO

@connect
@store
def storeTreeByTaxid(tree, taxid, *args, **kwargs):
    return  f"tree:{taxid}", tree.serialize()

@connect
@store
def storeVectorByTaxid(tree, taxid, *args, **kwargs):
    return  f"vector:{taxid}", tree.vectorize()

@connect
@store
def storeCulledVector(vector, taxid, cmin, cmax, fmax, *args, **kwargs):
    return  f"_vector:{taxid}:{cmin}:{cmax}:{fmax}", vector
    
@connect
@delete
def delTreeByTaxids(taxids, *args, **kwargs):
    return [f"tree:{_}" for _ in taxids]

@connect
@delete
def delVectorByTaxids(taxids, *args, **kwargs):
    return [f"vector:{_}" for _ in taxids]

@connect
@get
def getUniversalTree(taxid, *args, raw=False, **kwargs):
    return (f"tree:{taxid}", loadUnivGO)

@connect
@get
def getUniversalVector(taxid, *args, raw=False, **kwargs):
    return (f"vector:{taxid}", None)

@connect
@get
def getCulledVector(taxid, cmin, cmax, fmax, *args, raw=False, **kwargs):
    return  (f"_vector:{taxid}:{cmin}:{cmax}:{fmax}", None)

@connect
@listKey
def listTreeKey(*args, prefix=False, **kwargs):
    return ('tree:*', 'tree:')

@connect
@listKey
def listVectorKey(*args, prefix=False, **kwargs):
    return ('vector:*', 'vector:')

@connect
@listKey
def listCulledVectorKey(*args, prefix=False, **kwargs):
    return ('_vector:*', '_vector:')

def clearVectors(*args, **kwargs):
    keyList = [ _ for _ in listVectorKey() ]
    print(f"Clearing following vector elements content:\n\t{keyList}")
    delVectorByTaxids(keyList)
    
def clearTrees(*args, **kwargs):
    keyList = [ _ for _ in listTreeKey() ]
    print(f"Clearing following tree elements content:\n\t{keyList}")
    delTreeByTaxids(keyList)
    
def clear(*args, **kwargs):
    clearVectors(*args, **kwargs)
    clearTrees(*args, **kwargs)

"""
scan 0 MATCH *11*
"""

"""
def delTreeByTaxids(taxids):
    miss = []
    print(f"delete redis {taxids}")
    r = redis.Redis(host=HOST, port=PORT, db=0)       
    for taxid in taxids:
        _ = r.delete(f"tree:{taxid}")
        print(f"delete redis status {str(_)}")
        if int(str(_)) != 1:
            miss.append(taxid)
    return miss:
    
def storeTreeByTaxid(tree, taxid):
    print(f"redis taxid storage adding {taxid} {tree}")
    r = redis.Redis(host=HOST, port=PORT, db=0)
    if r.exists(f"tree:{taxid}"):
        raise KeyError(f"StoreTree error: taxid {taxid} already exists in store")
    d = tree.serialize()
    r.set(taxid, json.dumps(d))

def storeVectorByTaxid(vector, taxid):
    print(f"redis taxid storage adding {taxid} {vectors}")
    r = redis.Redis(host=HOST, port=PORT, db=0)
    if r.exists(f"vector:{taxid}"):
        raise KeyError(f"StoreVector error: taxid {taxid} already exists in store")
    
    r.set(taxid, json.dumps(vector))

def getUniversalTree(taxid, raw=False):
    print(f"redis taxid storage getting {taxid}")
    r = redis.Redis(host=HOST, port=PORT, db=0)
    _ = r.get(f"tree:{taxid}")
    if not _:
        raise KeyError(f"No taxid {taxid} found in tree store")
    # _ is bytes
    return _ \
        if raw \
        else loadUnivGO( json.loads(_) )

"""