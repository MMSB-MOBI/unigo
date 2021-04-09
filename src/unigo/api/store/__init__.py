from flask import Flask, jsonify, abort, request, make_response, Response
from ..data_objects import CulledGoParametersSchema as goParameterValidator
from ..data_objects import loadUnivGO
from .cache import setCacheType, delTreeByTaxids, storeTreeByTaxid, getTaxidKeys, getUniversalVector, getUniversalTree
from .cache import wipe as wipeStore
from .cache import unBuildtTreeIter, listTrees, listVectors
#GO_API_PORT = 1234
import time

from multiprocessing import Process, Value, Semaphore

#tasks = {
#    'build': Value('i', 0)
#}
bSemaphore = None
_MAIN_ = False

def bootstrap(trees=None, taxids=None, cacheType='local',\
    wipe=False, _main_=False, **kwargs):
    global _MAIN_
    print("Listening")
    print(trees)
    print(taxids)

    if _main_:
        print("Bootstraping main process")
        _MAIN_ = _main_
        global bSemaphore
        bSemaphore = Semaphore(1) # bleeding eyes

    setCacheType(cacheType, **kwargs)
    if wipe:
        wipeStore()


    for tree, taxid in zip(trees, taxids):
        storeTreeByTaxid(tree, taxid)

    app = Flask(__name__)

    app.add_url_rule('/handshake', 'handshake', handshake)
    
    app.add_url_rule('/taxids', 'view_taxids', view_taxids)
    
    app.add_url_rule('/unigo/<taxid>', 'view_unigo', view_unigo, methods=['GET'])
    
    app.add_url_rule('/vector/<taxid>', 'view_vector', view_vector)
    
    app.add_url_rule("/vector/<taxid>", 'view_culled_vector', view_culled_vector, methods=['POST'])
    
    app.add_url_rule('/add/unigo/<taxid>', 'add_unigo', add_unigo, methods=['POST'])

    app.add_url_rule('/del/unigo/<taxid>', 'del_unigo', del_unigo, methods=['DELETE'])

    app.add_url_rule('/build/vectors', 'build_vectors', build_vectors)

    app.add_url_rule('/list/<elemType>', 'list_elements', list_elements)

    return app

def list_elements(elemType):
    if elemType == 'vectors':
        return jsonify({ 'vectors': listVectors() })
      
    elif elemType == 'trees':
        return jsonify({ 'trees': listTrees() })

    print(f"Unknwon element type {elemType} to list")
    abort(404)

# vectorized all non-vectorized trees
def build_vectors():
    #print(f"bSemaphore value is {bSemaphore}")
    
   # if tasks['build'].value == 0:
   #     tasks['build'].value = 1
    if _MAIN_:
        # global bSemaphore exists
        if bSemaphore.acquire(block=False):
            print(f"bSemaphore is acquired")
            p = Process(
                target=my_func,#unBuildtTreeIter,
                args=(bSemaphore,),
                daemon=True
            )
            p.start()
            return jsonify({"status": "starting"}), 200
            
        print(f"bSemaphore is locked")  
            #tasks['build'].join()
            #tasks['build'] = None

        return jsonify({"status": "running"}), 202
    else: 
        print("##### Twilight zone ######") 
    #unBuildtTreeIter()
#sem = asyncio.Semaphore(10)


def my_func(_bSemaphore):
    time.sleep(10)
    print("Process finished")
    #num.value = 0
    _bSemaphore.release()
    print(f"Relasing bSemaphore")
def handshake():
    return "Hello world"

def view_taxids():
    return str( getTaxidKeys() )

def view_unigo(taxid):
    try:
        tree = getUniversalTree(taxid)
        return jsonify(tree.serialize())
    except KeyError as e:
        print(e)
        abort(404)

def view_vector(taxid):
    try:
        taxidVector = getUniversalVector(taxid)
    except KeyError as e:
        print(e)
        abort(404)
    
    return jsonify(taxidVector)

# Do we store culled vector ?
# vector:taxid:ccmin:cmax:fmax

def view_culled_vector(taxid):
    try:
        taxidVector = getUniversalVector(taxid)
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

    print(data)
    print(goParameter)
    print(taxidVector)
    try:
        _ = getCulledVector(taxid, cmin, cmax, fmax)
    except KeyError : # Culled vector not in cache, build&store

        _ = { # Do we got it all ? Dump it once
            'registry' : taxidVector['registry'],
            'terms' : {  goID : goTerm for goID, goTerm in taxidVector['terms'].items()\
                                                    if len(goTerm['elements']) >= cmin and\
                                                       len(goTerm['elements']) <= cmax and\
                                                       goTerm['freq']          <= fmax \
            }
        }
        print(f"look{_}")
        storeCulledVector(_, taxid, cmin, cmax, fmax)

    return jsonify(_)

def add_unigo(taxid):
    """ Add a unigo object through client API
        Parameters:
        ---------- 
        taxid: ncbi taxid 
        json payload: ..Univgo.serialize() dict
    """
    try:
        univGoTree = loadUnivGO(request.get_json())
        storeTreeByTaxid(univGoTree, taxid)
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
  
    return jsonify({"taxid" : "deletion OK"})