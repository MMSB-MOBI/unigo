from flask import Flask, jsonify, abort

UNIVERSAL_TREES = {}

UNIVERSAL_VECTORS = {}

def listen(trees=None, taxids=None):
    print("Listening")
    print(trees)
    print(taxids)
    global UNIVERSAL_TREES

    for txs, tr in zip(taxids, trees):
        for tx in txs:
            if tx in UNIVERSAL_TREES:
                raise KeyError(f"{tx} tree already exists")
        
            UNIVERSAL_TREES[tx] = tr
    
    print(UNIVERSAL_TREES)
    app = Flask(__name__)

    app.add_url_rule('/', 'index', index)
    
    app.add_url_rule('/taxids', 'view_taxids', view_taxids)
    
    app.add_url_rule('/unigo/<taxid>', 'view_unigo', view_unigo)
    
    app.add_url_rule('/vector/<taxid>', 'view_vector', view_vector)
    
    return app


def index():
    return "Hello world"

def view_taxids():
    global UNIVERSAL_TREES
    return str( list(UNIVERSAL_TREES.keys()) )


def view_unigo(taxid):
    print(f"/unigo/{taxid}")
    global UNIVERSAL_TREES
    if taxid in UNIVERSAL_TREES:
        d = UNIVERSAL_TREES[taxid].serialize()
        return jsonify(d)
    abort(404)

def view_vector(taxid):
    print(f"/vector/{taxid}")
    global UNIVERSAL_TREES
    global UNIVERSAL_VECTORS

    if taxid in UNIVERSAL_TREES:
        if taxid not in UNIVERSAL_VECTORS:
            print(f"Creating vectors for {taxid}")
            UNIVERSAL_VECTORS[taxid] = UNIVERSAL_TREES[taxid].vectorize()
        else:
            print(f"Found in vector cache {taxid}")

        return jsonify(UNIVERSAL_VECTORS[taxid])
    abort(404)

