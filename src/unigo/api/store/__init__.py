from flask import Flask, jsonify, abort, request
from ..data_objects import CulledGoParametersSchema as goParameterValidator
from ..data_objects import loadUnivGO
from .cache import setCacheType, delTreeByTaxids, storeTreeByTaxid, getTaxidKeys, getUniversalVector, getUniversalTree
from .cache import wipe as wipeStore
#GO_API_PORT = 1234

def bootstrap(trees=None, taxids=None, cacheType='local',\
    wipe=False, **kwargs):
    print("Listening")
    print(trees)
    print(taxids)
   
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

    return app

# vectorized all non-vectorized trees
def build():
    pass

def handshake():
    return "Hello world"

def view_taxids():
    return str( getTaxidKeys() )


def view_unigo(taxid):
    print(f"/unigo/{taxid}")
    try:
        tree = getUniversalTree(taxid)
        return jsonify(tree.serialize())
    except KeyError as e:
        print(e)
        abort(404)

        

def view_vector(taxid):
    print(f"/vector/{taxid}")

    try:
        taxidVector = getUniversalVector(taxid)
    except KeyError as e:
        print(e)
        abort(404)
    
    return jsonify(taxidVector)
    
def view_culled_vector(taxid):
    try:
        taxidVector = getUniversalVector(taxid)
    except KeyError as e:
        print(e)
        abort(404)

    data = request.get_json()
    _goParameterValidator = goParameterValidator()
    goParameter = _goParameterValidator.load(request)
    print(data)
    print(goParameter)
    print(taxidVector)

    _ = {
        'registry' : taxidVector['registry'],
        'terms' : {  goID : goTerm for goID, goTerm in taxidVector['terms'].items()\
                                                    if len(goTerm['elements']) >= goParameter['minCount'] and\
                                                       len(goTerm['elements']) <= goParameter['maxCount'] and\
                                                       goTerm['freq']          <= goParameter['maxFreq']\
        }
    }

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
        print(f"{taxid} already exist in database, reject insertion")
        abort(403)
    except Exception as e:
        print(e)
        abort(500)

def del_unigo(taxid):
    print(f"del_unigo {taxid}")

    try:
        delTreeByTaxids([taxid])
    except KeyError:
        print(f"{taxid} not found in database, nothing to delete")       
        abort(404) 
  
    return jsonify({"taxid" : "deletion OK"})