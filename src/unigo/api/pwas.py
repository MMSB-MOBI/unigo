from flask import Flask, abort, jsonify, request
import requests, json
from marshmallow import EXCLUDE
from .. import vloads as createGOTreeTestUniverseFromAPI
from .. import uloads as createGOTreeTestFromAPI
from .. import utils
from .io import checkPwasInput
from ..stat_utils import applyOraToVector, kappaClustering
from .data_objects import CulledGoParametersSchema as goParameterValidator
from copy import deepcopy as copy

GOPORT=1234
GOHOST="127.0.0.1"
def listen(goApiHost:str, goApiPort:int, vectorized:bool):
    global GOPORT, GOHOST
    GOPORT = goApiPort
    GOHOST = goApiHost

    app = Flask(__name__)
    app.config['JSON_SORT_KEYS'] = False # To keep dict order in json
    app.add_url_rule("/", 'hello', hello)
    if vectorized:
        print(f"PWAS API vector listening")
        app.add_url_rule("/compute", "computeOverVector", computeOverVector, methods=["POST"])
    else:
        print(f"PWAS API tree listening")
        app.add_url_rule("/compute", "computeOverTree", computeOverTree, methods=["POST"])
    
    app.add_url_rule("/loadVector/<taxid>", loadVector)
    return app

def hello():
    return "Hello pwas"

def computeOverVector():
    forceUniversal = False
    data = checkPwasInput()
   
    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')
    
    if forceUniversal:
        go_resp = utils.unigo_vector_from_api(GOHOST, GOPORT, data["taxid"])
    else:
         # Culling vector parameters
        _goParameterValidator = goParameterValidator()
        goParameter = _goParameterValidator.load(request,  unknown=EXCLUDE)
        go_resp = utils.unigo_culled_from_api(GOHOST, GOPORT, data["taxid"], goParameter)

    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    
    # Here Contract of 3 NS on go_resp
    vectorizedProteomeTrees = go_resp.json()
 #  print(vectorizedProteomeTrees.keys())
    
    kappaClusters = kappaClusteringOverNS(vectorizedProteomeTrees, data, merge=True)
    return jsonify(kappaClusters)

def fuseVectorNameSpace(_vectorElements, merge):
    vectorElements = copy(_vectorElements)

    if not merge:
        for ns, vectorElement in vectorElements.items():
            for goID, goVal in vectorElement["terms"].items():
                goVal["ns"] = ns
        return vectorElements

    fusedNS = {
        "terms" : {},
        "registry": []
    }

    for ns, vectorElement in vectorElements.items():
        # All registries are identical
        if not fusedNS["registry"]:
            fusedNS["registry"] = vectorElement["registry"]
        for goID, goVal in vectorElement["terms"].items():
            goVal["ns"] = ns
            fusedNS["terms"][goID]  = goVal
    return { "fusedNS" : fusedNS }


def kappaClusteringOverNS(_vectorElements, expData, merge=False):

    vectorElements = fuseVectorNameSpace(_vectorElements, merge)  
    kappaClusters = {}   
    for ns, vectorElement in vectorElements.items():
        res = applyOraToVector(vectorElement,\
            expData["all_accessions"],\
            expData["significative_accessions"],\
            0.05)
        if len( res.keys() ) <= 1:
            kappaClusters[ns] = res
        else:
            Z = kappaClustering(vectorElement["registry"], res)
            kappaClusters[ns] = Z
    
    return(kappaClusters)

def computeOverTree():
    data = checkPwasInput() 
    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')

    go_resp = utils.unigo_tree_from_api(GOHOST, OPORT, data["taxid"])

    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    
    print("Create Unigo tree")
    unigoTreeFromAPI = createGOTreeTestFromAPI(go_resp.text, data["all_accessions"])
    x,y = unigoTreeFromAPI.dimensions
    print("Unigo Object successfully buildt w/ following dimensions:")
    print(f"\txpTree => nodes:{x[0]} children_links:{x[1]}, total_protein_occurences:{x[2]}, protein_set:{x[3]}")  
    print(f"\t universeTree => nodes:{y[0]} children_links:{y[1]}, total_protein_occurences:{y[2]}, protein_set:{y[3]}")  

    #Compute fisher stats
    if data["method"] == "fisher":
        print("Computing ORA with fisher")
        rankingsORA = unigoTreeFromAPI.computeORA(data["significative_accessions"], verbose = False)
        
        return rankingsORA.json

    return {"not computed": "unavailable stat method"}

def _loadVector(taxid):
    go_resp = utils.unigo_vector_from_api(GOPORT, taxid)
    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    else:
        return go_resp

def loadVector(taxid):
    if _loadVector(taxid):
        return {"ok" : f"Vector loaded for {taxid} taxid"}
    else:
        return {"error" : f"Can't load vector for {taxid} taxid"}




