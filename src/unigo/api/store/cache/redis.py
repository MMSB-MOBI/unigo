from pyrediscore import connect, store, get, delete, listKey, setDatabaseParameters
from ...data_objects import loadUnivGO
from ....tree import enumNSkeys, assertAndCoherceValidNamespace
from decorator import decorator


######## SETTER ##########
@connect
@store
def storeTreeByTaxid(tree, taxid, *args, **kwargs):
    ns =  enumNSkeys[tree.ns]
    return  f"tree:{taxid}:{ns}", tree.serialize()

@connect
@store
def storeVector(tree, taxid=None, key=None, *args, **kwargs):
    if key:
        return  f"vector:{key}", tree.vectorize()
    if taxid:
        ns =  enumNSkeys[tree.ns]
    return  f"vector:{taxid}:{ns}", tree.vectorize()

    raise ValueError(f"No Key can be deduced for {tree}")

def storeVectorByTaxid(tree, taxid, *args, **kwargs):
    _ = storeVector(tree, taxid=taxid, *args, **kwargs)
    return _

@connect
@store
def storeCulledVector(vector, taxid, ns, cmin, cmax, fmax, *args, **kwargs):
    #print(f"##{ns}")
    #ns =  enumNSkeys[ vector["ns"] ]
    return  f"_vector:{taxid}:{ns}:{cmin}:{cmax}:{fmax}", vector


@decorator
def setNamespace(fn, *args, ns=None, **kwargs):
    _ns = list(enumNSkeys.values())
    if ns:
        _ = assertAndCoherceValidNamespace(ns)
        _ns = [_]
    return fn(*args, ns=_ns, **kwargs)

######## DELETE ##########
# These operations will complete but may raise key error if at least one entry is missing
@connect
@delete
@setNamespace
def delTreeByTaxids(taxids, *args, ns=None, **kwargs):
    """ Delete the 3 NS universal trees of the provided taxids 
    option: ns parameter can be used to delete a single namespace tree
    """
    #print("BIPBIP") 
    #print(ns)
    #print(taxids)
    print("==>", [f"tree:{taxid}:{currNS}" for currNS in ns for taxid in taxids] )
    
    return [f"tree:{taxid}:{currNS}" for currNS in ns for taxid in taxids]

@connect
@delete
@setNamespace
def delVectorByTaxids(taxids, *args, ns=None, **kwargs):
    """ Delete the 3 NS universal vectors of the provided taxids 
    option: ns parameter can be used to delete a single namespace tree
    """

    return [f"vector:{taxid}:{currNS}" for currNS in ns for taxid in taxids]

@connect
@delete # No @setNameSpace as enumeration of culled requires list operation
def delCulledVectorByTaxids(taxids, *args, ns=None, **kwargs):
    """ Delete the all culled vectors under provided taxids 
    option: ns parameter can be used to delete the set of culled vectors under
    this namespace
    """
    l = []
    for t in taxids:
        l += listCulledVectorKey(taxid=t, ns=ns)
    return [f"_vector:{_}" for _ in l]

@connect
@delete
def delVectorByKeys(keys, *args, **kwargs):
    return [f"vector:{key}" for key in keys]

@connect
@delete
def delTreeByKeys(keys, *args, **kwargs):
    return [f"tree:{key}" for key in keys]

@connect
@delete
def delCulledVectorByKeys(keys, *args, **kwargs):
    return [f"_vector:{key}" for key in keys]

def deleteTaxids(taxids, *args, **kwargs):
    deletedT = delTreeByTaxids(taxids, ignore=True)
    deletedV = delVectorByTaxids(taxids, ignore=True)
    deletedC = delCulledVectorByTaxids(taxids, ignore=True)
    print('TC:: ', deletedT, deletedV, deletedC)
    if deletedC == 0 and deletedV == 0 and deletedT == 0:
        return None
    return (deletedT, deletedV, deletedC)
######## GETTER ##########

@setNamespace
def getUniversalTree3NS(taxid, *args, ns=None, raw=False, **kwargs):
    return { _ns : getUniversalTree(f"{taxid}:{_ns}", *args, raw=raw, **kwargs)\
              for _ns in ns }
        
@connect
@get
def getUniversalTree(fullKey, *args, raw=False, **kwargs):   
    return (f"tree:{fullKey}", loadUnivGO)
   
@setNamespace
def getUniversalVector3NS(taxid, *args, ns=None, raw=False, **kwargs):
    return { _ns: getUniversalVector(taxid, _ns, *args, raw=raw, **kwargs)\
              for _ns in ns }
    
@connect
@get
def getUniversalVector(taxid, ns, *args, raw=False, **kwargs):
    return (f"vector:{taxid}:{ns}", None)

@setNamespace
def getCulledVector3NS(taxid, cmin, cmax, fmax, *args, ns=None, raw=False, **kwargs):
    return {_ns : getCulledVector(taxid, _ns, cmin, cmax, fmax, *args, raw=False, **kwargs)\
              for _ns in ns }
    #return  (f"_vector:{taxid}:{cmin}:{cmax}:{fmax}", None)

@connect
@get
def getCulledVector(taxid, ns, cmin, cmax, fmax, *args, raw=False, **kwargs):    
    return (f"_vector:{taxid}:{ns}:{cmin}:{cmax}:{fmax}", None)


######## LIST ##########
@decorator
def expandQueryList(fn, *args, **kwargs):
    #print(kwargs)
    query = [ kwargs['taxid'] if not kwargs['taxid'] is None else  '*',\
              kwargs['ns']    if not kwargs['ns'] is None else  '*'
    ]
    if 'cmin' in kwargs: # One test to assume it is a culled related search
        query += [
              kwargs['cmin']  if not kwargs['cmin'] is None else  '*',\
              kwargs['cmax']  if not kwargs['cmax'] is None else  '*',\
              kwargs['fmax']  if not kwargs['fmax'] is None else  '*',\
        ]
    
    _ = ':'.join(query) 
    return fn(_, *args, **kwargs)

@connect
@listKey
@expandQueryList
def listTreeKey(*args, taxid=None, ns=None, prefix=False, **kwargs):
    #return ('tree:*', 'tree:')
    #print(f"LIST TREE KEY PREFIX {args[0]}")
    return (f"tree:{args[0]}", 'tree:')

@connect
@listKey
@expandQueryList
def listVectorKey(*args, taxid=None, ns=None, prefix=False, **kwargs):
    #return ('vector:*', 'vector:')
    #print(f"LIST VECTOR KEY PREFIX {args[0]}")
    return (f"vector:{args[0]}", 'vector:')

@connect
@listKey
@expandQueryList
def listCulledVectorKey(*args, taxid=None, ns=None, cmin=None, cmax=None, fmax=None, prefix=False, **kwargs):
    #print(f"LIST CULLED KEY PREFIX {args[0]}")
    #return ('_vector:*', '_vector:')
    return(f"_vector:{args[0]}", '_vector:')

######### CLEANUP #########
# LIST base cleanup, deletion operation can't fail
def clearTrees(*args, **kwargs):
    keyList = [ _ for _ in listTreeKey() ]
    print(f"Clearing following tree elements content:\n\t{keyList}")
    delTreeByKeys(keyList)
    return len(keyList)

def clearVectors(*args, **kwargs):
    keyList = [ _ for _ in listVectorKey() ]
    print(f"Clearing following vector elements content:\n\t{keyList}")     
    delVectorByKeys(keyList)
    return len(keyList)

def clearCulledVectors(*args, **kwargs):
    keyList = [ _ for _ in listCulledVectorKey() ]
    print(f"Clearing following culled vectors elements content:\n\t{keyList}")
    delCulledVectorByKeys(keyList)
    return len(keyList)
    
def clear(*args, **kwargs):
    print("Clearing Entiere database content")
    _ = clearTrees(*args, **kwargs)
    print(f"{_} tree(s) cleared")
    _ = clearVectors(*args, **kwargs)
    print(f"{_} vector(s) cleared")
    _ = clearCulledVectors(*args, **kwargs)
    print(f"{_} culled vector(s) cleared")

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