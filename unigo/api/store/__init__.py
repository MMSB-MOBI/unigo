from flask import Flask, jsonify, abort, request, make_response, Response
from ..data_objects import CulledGoParametersSchema as goParameterValidator
from ... import Unigo
from .cache import setCacheType, delTreeByTaxids, storeTreeByTaxid, getTaxidKeys, getUniversalTrees
from .cache import getCulledVectors, storeCulledVector, getUniversalVectors, delVectorsByTaxid
from .cache import clear as clearStore, deleteTaxids
from .cache import buildUniversalVector, listTrees, listVectors, listMissUniversalVector, listCulled
from .cache import status as storeStatus
from decorator import decorator
from ...tree import enumNSkeysRevert as humanReadableNS
# IF OFFLOADABLE BLOCKING ROUTES KEEP ON SHOWING UP
# https://blog.miguelgrinberg.com/post/using-celery-with-flask

import time

from multiprocessing import Process, Value, Semaphore

bSemaphore = None
C_TYPE = None
_MAIN_ = False

def bootstrap(newElem=None, cacheType='redis',\
    clear=False, _main_=False, **kwargs):
    global _MAIN_, C_TYPE
    
    if _main_:
        print("Bootstraping main process")
        _MAIN_ = _main_
        global bSemaphore
        bSemaphore = Semaphore(1) # bleeding eyes
        setCacheType(cacheType, **kwargs)
        C_TYPE = cacheType
        if clear:
            clearStore()

        if newElem:
            print("Boostraping unigo store with following new ressources")
            for taxid, ns, tree in newElem:
                print(f"\tAdding{taxid}:{tree}")
                storeTreeByTaxid(tree, taxid)

        tcount, vcount, _vcount = storeStatus()
        print(f"Database content:\n\t{tcount} trees, {vcount} vectors, {_vcount} culled")

        app = Flask(__name__)

        app.add_url_rule('/ping', 'ping', ping)
        
        # TO FIX
        #app.add_url_rule('/taxids', 'view_taxids', view_taxids)
        
        app.add_url_rule('/unigos/<taxid>', 'view_unigos', view_unigos, methods=['GET'])
        
        app.add_url_rule('/vectors/<taxid>', 'view_vectors', view_vectors, methods=['GET'])
        
        app.add_url_rule("/vectors/<taxid>", 'view_culled_vectors', view_culled_vectors, methods=['POST'])
        
        app.add_url_rule('/add/taxid/<taxid>', 'add_unigo3NS', add_unigo3NS, methods=['POST'])

        app.add_url_rule('/del/taxid/<taxid>', 'del_taxonymy', del_taxonymy, methods=['DELETE'])

        app.add_url_rule('/add/unigo/<taxid>', 'add_unigo', add_unigo, methods=['POST'])

        app.add_url_rule('/del/unigo/<taxid>', 'del_unigo', del_unigo, methods=['DELETE'])

        app.add_url_rule('/build/vectors', 'build_vectors', build_vectors)

        app.add_url_rule('/list/<elemType>', 'list_elements', list_elements)

        return app

def list_elements(elemType):
    if elemType == 'vectors':
        _ = listVectors()
        print(_)
        return jsonify({ 'vectors':  _ })
      
    elif elemType == 'trees':
        return jsonify({ 'trees': listTrees() })
    
    elif elemType == 'culled':
        return jsonify({ 'culled': listCulled() })


    print(f"Unknwon element type {elemType} to list")
    abort(404)

# CPU intensive route -> offloading it to multiprocess
# using semaphore to prevent concurrency popup on multiple requests
def build_vectors():
    if _MAIN_:
        missUniversalVector = listMissUniversalVector()
        if not missUniversalVector:
            return jsonify({"status": "nothing to build"}), 200

        # global bSemaphore exists
        if bSemaphore.acquire(block=False):
            print(f"bSemaphore is acquired")
            p = Process(
                target=_buildUniversalVector,
                args=(bSemaphore, C_TYPE),
                daemon=True
            )
            p.start()
            return jsonify({
                "status": "starting",
                "targets" : missUniversalVector
            }), 200
        
        print(f"bSemaphore is locked")   
        return jsonify({"status": "running", "targets" : missUniversalVector}), 202
    else: 
        print("##### Twilight zone ######") 


@decorator
def semHolder(fn, _bSemaphore, cacheType, *args, **kwargs):
    #print("Decrorator start")
    setCacheType(cacheType)
    #time.sleep(10)
    fn(*args, **kwargs)
    _bSemaphore.release()
    #print(f"Releasing bSemaphore")

@decorator
def nsHumanizer(fn, *args, **kwargs):
    """ request.json wraps dict index by first translating its NS one-key in long description """
    _ = fn(*args, **kwargs)
    return jsonify( {
            humanReadableNS[k] : v for k, v in _.items()            
        } )

@semHolder
def _buildUniversalVector(*args,**kwargs):
    buildUniversalVector()

def ping():
    return "pong"

def view_taxids():
    return str( getTaxidKeys() )

@nsHumanizer
def view_unigos(taxid, *args, **kwargs):
    try:
        trees = getUniversalTrees(taxid)
    except KeyError as e:
        print(f"Key Error at unigo identifer (aka: taxid) \"{taxid}\"\n=>{e}")
        abort(404)
    return { ns : tree.serialize() for ns, tree in trees.items() }
##################################################
# ACCESSORS to taxid specific goterms as vectors
##################################################
# GET route, return universal vector # ADD NS in route ?
#@nsHumanizer
@nsHumanizer
def view_vectors(taxid, *args, **kwargs):
    try:
        taxidVectors = getUniversalVectors(taxid)
    except KeyError as e:
        print(e)
        abort(404)
    
    return taxidVectors

# POST route
# Presuming live culling is not CPU intensive enough to
# promote multiprocess offloading:: CHECK IT
@nsHumanizer
def view_culled_vectors(taxid, *args, **kwargs): 
    try:
        taxidVectors = getUniversalVectors(taxid)
    except KeyError as e:
        print(e)
        abort(404)
    try:
        data = request.get_json()
        _goParameterValidator = goParameterValidator()
        goParameter = _goParameterValidator.load(request)
        cmin, cmax, fmax = ( goParameter['minCount'],\
            goParameter['maxCount'], goParameter['maxFreq'])
    except Exception as e:
        print(e)
        print(f"Malformed GO parameters\n=>{data}")
        abort(400)
    print("Searching for culled vector...")
    try:
        culledVectors = getCulledVectors(taxid, cmin, cmax, fmax)
        print(f"Culled vector triplet found ! [{taxid} {cmin} {cmax} {fmax}]")
    except KeyError : # Culled vector not in cache, build&store
        print(f"Building culled vectors... [{taxid} {cmin} {cmax} {fmax}]")
        culledVectors = {}
        for tree_ns_key, taxidVector in taxidVectors.items():
            _ = { # Do we got it all ? Dump it once
                'registry' : taxidVector['registry'],
                'terms' : {  goID : goTerm for goID, goTerm in taxidVector['terms'].items()\
                                                        if len(goTerm['elements']) >= cmin and\
                                                           len(goTerm['elements']) <= cmax and\
                                                           goTerm['freq']          <= fmax \
                }
            }            
            storeCulledVector(_, taxid, tree_ns_key, cmin, cmax, fmax)
            culledVectors[tree_ns_key] = _

    return culledVectors



def add_unigo3NS(taxid):
    """ Add 3NS unigo object through client API
        Parameters:
        ---------- 
        taxid: ncbi taxid 
        json payload: { NS : Univgo.serialize(), }
    """
    print(f"Adding unigo3NS for taxid: {taxid}")
    _ = request.get_json()
    for tree in _.values():
        try:
            unigo_blueprint = Unigo(from_serial=tree)
           
            storeTreeByTaxid(unigo_blueprint, taxid)            
        except KeyError as e:
            print(f"Similar tree key already exist in database or malformed request, reject insertion\n=>{e}")
            abort(403)
        except Exception as e:
            print(f"add unigo internal error:{e}")
            abort(500)
    return jsonify({ taxid + " 3NS" : "insertion OK"})

def del_taxonymy(taxid):
    _ = deleteTaxids([taxid]) #delTreeByTaxids([taxid])
    if _ is None: 
        print(f"{taxid} not found in database, nothing to delete")       
        abort(404) 
    
    return jsonify({f"{taxid}" : { "deletions [Tree, Vector, Culled]" : list(_) }             
            })

# ADD Subroutes to delete culled and/or vector whilre preserving Trees ?

## MAY BE DECLASSIFIED
def add_unigo(taxid):
    """ Add a unigo object through client API
        Parameters:
        ---------- 
        taxid: ncbi taxid 
        json payload: ..Univgo.serialize() dict
    """
    try:
        unigo_blueprint= Unigo(from_serial = request.get_json())
        storeTreeByTaxid(unigo_blueprint, taxid)
        return jsonify({"taxid" : "insertion OK"})
    except KeyError as e:
        print(f"{taxid} already exist in database or malformed request, reject insertion\n=>{e}")
        abort(403)
    except Exception as e:
        print(f"add unigo internal error:{e}")
        abort(500)

def del_unigo(taxid):
    try:
        delTreeByTaxids([taxid])
    except KeyError:
        print(f"{taxid} not found in database, nothing to delete")       
        abort(404) 
    
    print(f"deleting related vectors...")
    delVectorsByTaxid([taxid])

    return jsonify({"taxid" : "deletion OK"})

