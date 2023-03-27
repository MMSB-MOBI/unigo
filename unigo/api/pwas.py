from flask import Flask, abort, jsonify, request
from marshmallow import EXCLUDE
from .. import Unigo, unigo_obs_mask

from .store.client.viewers import unigo_tree_from_api, unigo_vector_from_api, unigo_culled_from_api
from .io import checkPwasInput
from ..stats.clustering import kappaClustering
from ..stats.ora import apply_ora_to_unigo, applyOraToVector
from .data_objects import CulledGoParametersSchema as goParameterValidator
from copy import deepcopy as copy
from .store.client import handshake

def listen(vectorized:bool):
    
    app = Flask(__name__)
    app.config['JSON_SORT_KEYS'] = False # To keep dict order in json
    app.add_url_rule("/ping", 'hello', hello)
    if vectorized:
        print(f"PWAS API vector listening")
        app.add_url_rule("/compute", "computeOverVector", computeOverVector, methods=["POST"])
    else:
        print(f"PWAS API tree listening")
        app.add_url_rule("/compute", "computeOverTree", computeOverTree, methods=["POST"])
    
    app.add_url_rule("/loadVector/<taxid>", loadVector) # WARNING : don't work
    app.add_url_rule('/get_members', 'get_members', get_members, methods=['POST'])
    return app

def hello():
    return "Hello pwas"

def computeOverVector():
    forceUniversal = False
    data = checkPwasInput()
   
    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')
    
    if forceUniversal:
        go_resp = unigo_vector_from_api(data["taxid"])
    else:
         # Culling vector parameters
        _goParameterValidator = goParameterValidator()
        goParameter = _goParameterValidator.load(request,  unknown=EXCLUDE)
        #go_resp = unigo_culled_from_api(GOHOST, GOPORT, data["taxid"], goParameter)
        go_resp = unigo_culled_from_api(data["taxid"], goParameter)
    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    
    # Here Contract of 3 NS on go_resp
    vectorizedProteomeTrees = go_resp.json()

    kappaClusters = kappaClusteringOverNS(vectorizedProteomeTrees, data)
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

    #print("expData", expData)
    if expData.get("pvalue"):
        pvalue = expData["pvalue"]
    else:
        pvalue = 0.05

    print("pvalue", pvalue)

    for ns, vectorElement in vectorElements.items():
        res = applyOraToVector(vectorElement,\
            expData["all_accessions"],\
            expData["significative_accessions"],\
            pvalue)
        formatted_res = [{**{"go": go_term}, **res[go_term]} for go_term in res]
        if len( res.keys() ) <= 1:
            kappaClusters[ns] = {'Z': res, 'list' : formatted_res}
        else:
            Z = kappaClustering(vectorElement["registry"], res)
            kappaClusters[ns] = {'Z': Z, 'list' : formatted_res}
    
    return(kappaClusters)

def computeOverTree():
    data = checkPwasInput() 
    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')

    _goParameterValidator = goParameterValidator()
    goParameter = _goParameterValidator.load(request,  unknown=EXCLUDE)
    #print(f"Incoming request GO constraints {goParameter}")
    #go_resp = unigo_tree_from_api(GOHOST, GOPORT, data["taxid"])
    go_resp = unigo_tree_from_api(data["taxid"])

    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    
    tree_elements = go_resp.json()
    ora_data_3NS = {}
    for ns, tree_element in tree_elements.items():
        blueprint_unigo = Unigo(from_serial=tree_element)
        mask_unigo      = unigo_obs_mask(blueprint_unigo, data["all_accessions"])
        dim = mask_unigo.dimensions
        print(f"Unigo Object \"{ns}\" successfully buildt and masked w/ following dimensions:")
        print(f"\t=> nodes:{dim[0]} children_links:{dim[1]}, total_protein_occurences:{dim[2]}, protein_set:{dim[3]}")  
        
    
    #Compute fisher stats
        if data["method"] == "fisher":
            print("Computing ORA with fisher")
            rankingsORA = apply_ora_to_unigo(mask_unigo,\
                            data["significative_accessions"],\
                            pvalue_max = 0.1, \
                            verbose = False)
            print(f"Filtering {ns} ORA analysis based on following GO constraints:{goParameter}")
            ora_data_3NS[ns] = [ go_ora for go_ora in rankingsORA\
                if go_ora["K_k_N_n_deep"][0] > goParameter['minCount'] and \
                   go_ora["K_k_N_n_deep"][0] <= goParameter['maxCount'] and \
                   go_ora["pathway_freq_deep"] <=   goParameter['maxFreq'] ]
    
    if data["method"] == "fisher":
        return ora_data_3NS
    return {"not computed": "unavailable stat method"}

def computeORA(unigo_tree, significative_accessions, pvalue, min_count, max_count, max_freq):
    rankingsORA = apply_ora_to_unigo(unigo_tree,\
                            significative_accessions,\
                            pvalue_max = pvalue, \
                            verbose = False)
    ora_data_3NS = {}
    for ns in unigo_tree:
        ora_data_3NS[ns] = [ go_ora for go_ora in rankingsORA\
            if go_ora["K_k_N_n_deep"][0] > min_count and \
                go_ora["K_k_N_n_deep"][0] <= max_count and \
                go_ora["pathway_freq_deep"] <=   max_freq ]
    
    return ora_data_3NS

def _loadVector(taxid):
    #go_resp = unigo_vector_from_api(GOPORT, taxid)
    go_resp = unigo_vector_from_api(taxid)
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
    
def get_members():
    data = check_get_members_input() 
    print("get_members route", data)
    go_resp = unigo_tree_from_api(data["collection"])
    print('go_resp', go_resp)
    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)
    
    tree_elements = go_resp.json()
    blueprint = Unigo(from_serial=tree_elements[data['ns']])
    
    response = {}
    for go in data['go_members']:
        response[go] = blueprint.getByID(go).getMembers()

    return response

def check_get_members_input():
    print('check get_members input')
    data = request.get_json()
    if not data:
        print(f"ERROR in request : empty posted data. Abort 400")
        abort(400)

    for key in ['collection', 'go_members', 'ns']:
        if not key in data:
            print(f"ERROR in request : {key} not in posted data. Abort 400")
            abort(400)

    return data




