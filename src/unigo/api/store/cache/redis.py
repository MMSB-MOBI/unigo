import redis
import json
from ...data_objects import loadUnivGO

HOST = 'localhost'
PORT = 6379

def wipe():
    print(f"Wiping redis stores content")
    r = redis.Redis(host=HOST, port=PORT, db=0)
    r.flushdb()

def delTreeByTaxids(taxids):
    miss = []
    print(f"delete redis {taxids}")
    r = redis.Redis(host=HOST, port=PORT, db=0)
    for taxid in taxids:
        _ = r.delete(taxid)
        print(f"delete redis status {str(_)}")
        if int(str(_)) != 1:
            miss.append(taxid)
    if miss:
        raise KeyError(f"{miss} to delete elements not found in redis store")

def setDatabaseParameters(host=HOST, port=PORT):
    global HOST, PORT
    HOST = host
    PORT = port
    print(f"Set Database parameters to {HOST}:{PORT}")

def storeTreeByTaxid(tree, taxid):
    print(f"redis taxid storage adding {taxid}")
    r = redis.Redis(host=HOST, port=PORT, db=0)
    if r.exists(taxid):
        raise KeyError(f"StoreTree error: taxid {taxid} already exists in store")
    d = tree.serialize()
    r.set(taxid, json.dumps(d))

def getUniversalTree(taxid, raw=False):
    print(f"redis taxid storage getting {taxid}")
    r = redis.Redis(host=HOST, port=PORT, db=0)
    _ = r.get(taxid)
    if not _:
        raise KeyError(f"No taxid {taxid} found in tree store")
    # _ is bytes
    return _ \
        if raw \
        else loadUnivGO( json.loads(_) )