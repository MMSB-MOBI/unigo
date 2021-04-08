import redis
import json, sys
from ...data_objects import loadUnivGO

from decorator import decorator

HOST = 'localhost'
PORT = 6379

def setDatabaseParameters(host=HOST, port=PORT):
    global HOST, PORT
    HOST = host
    PORT = port
    print(f"Set Database parameters to {HOST}:{PORT}")

def wipe():
    print(f"Wiping redis stores content")
    r = redis.Redis(host=HOST, port=PORT, db=0)
    r.flushdb()

@decorator
def connect(fn, *args, **kwargs):
    r = redis.Redis(host=HOST, port=PORT, db=0)
    return fn(r, *args, **kwargs)

@decorator
def delete(fn, r, *args, **kwargs):
    miss = []
    keys = fn(*args, **kwargs)
    for key in keys:
        _ = r.delete(key)
            #print(f"delete redis status {str(_)}")
        if int(str(_)) != 1:
            miss.append(taxid)
    if miss:
        raise KeyError(f"{miss} to delete elements not found in redis store")

@decorator
def store(fn, r, *args, **kwargs):
    print(args)
    key, obj = fn(*args, **kwargs)
    if r.exists(key):
        raise KeyError(f"Store error: {key} already exists in store")
    r.set(key, json.dumps(obj))

@decorator
def get(fn, r, *args, **kwargs):
    key, _deserializer = fn(*args, **kwargs)
    _ = r.get(key)

    if not _:
        raise KeyError(f"No key {key} found in store")

    if kwargs['raw']:
        return _
    d = json.loads(_)
    return d \
        if _deserializer is None \
        else _deserializer(d)

@decorator
def listKey(fn, r, *args, **kwargs):
    _pattern, _prefix = fn(*args, **kwargs)
    for hitKey in r.scan_iter(match=_pattern, count=None, _type=None):
        hitKey = hitKey.decode()
        yield hitKey if kwargs['prefix'] else hitKey.replace(_prefix, '')

@connect
@store
def storeTreeByTaxid(tree, taxid, *args, **kwargs):
    return  f"tree:{taxid}", tree.serialize()

@connect
@store
def storeVectorByTaxid(tree, taxid, *args, **kwargs):
    return  f"vector:{taxid}", tree.vectorize()
    
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
def getUniversalVector(taxid, raw=False):
    return (f"vector:{taxid}", None)

@connect
@listKey
def listTreeKey(*args, prefix=False, **kwargs):
    return ('tree:*', 'tree:')

@connect
@listKey
def listVectorKey(*args, prefix=False, **kwargs):
    return ('vector:*', 'vector:')

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