from .io import PwasData
from ..utils import unigo_vector_from_api, unigo_culled_from_api
from .data_objects import CulledGoParametersSchema as goParameterValidator
from marshmallow import EXCLUDE
from flask import abort, jsonify
from ..stat_utils import applyOraToVector, kappaClustering
from copy import deepcopy as copy
import numpy as np
import networkx as nx
from scipy.cluster.hierarchy import ward
from ..cluster_tree import create_go_clusters_tree


def compute_over_vector(data: PwasData, go_host:str, go_port:int, compute_clusters: bool = False, force_universal = True):
    print(f'I get data with {len(data["all_accessions"])} proteins accessions including {len(data["significative_accessions"])} significatives')

    if force_universal:
        go_resp = unigo_vector_from_api(go_host, go_port, data["taxid"])
    else:
         # Culling vector parameters
        _goParameterValidator = goParameterValidator()
        goParameter = _goParameterValidator.load(data,  unknown=EXCLUDE)
        go_resp = unigo_culled_from_api(go_host, go_port, data["taxid"], goParameter)

    if go_resp.status_code != 200:
        print(f"ERROR request returned {go_resp.status_code}")
        abort(go_resp.status_code)

    vectorizedProteomeTrees = go_resp.json()
    kappaClusters = pwas_compute_and_cluster(vectorizedProteomeTrees, data, compute_clusters)
    return kappaClusters


def pwas_compute_and_cluster(_vectorElements, expData, compute_clusters, merge=False):

    vectorElements = fuseVectorNameSpace(_vectorElements, merge)  
    kappaClusters = {}   

    #print("expData", expData)
    if expData.get("pvalue"):
        pvalue = expData["pvalue"]
    else:
        pvalue = 0.05

    print("pvalue", pvalue)

    for ns, vectorElement in vectorElements.items():
        print('compute', ns)
        res = applyOraToVector(vectorElement,\
            expData["all_accessions"],\
            expData["significative_accessions"],\
            pvalue)
        formatted_res = [{**{"go": go_term}, **res[go_term]} for go_term in res]
        kappaClusters[ns] = {'list' : formatted_res}
        if compute_clusters:
            threshold = 0.8
            GO_clusters = create_go_clusters_tree(formatted_res)
            kappaClusters[ns]['go_clusters'] = GO_clusters.serialize_to_d3()
            #kappaClusters[ns]["adjacency_matrix"] = {'threshold' : threshold, "matrix": adjacency_matrix}
            '''if len( res.keys() ) <= 1:
                kappaClusters[ns] = {'Z': res, 'list' : formatted_res}
            else:
                Z = kappaClustering(vectorElement["registry"], res)
                kappaClusters[ns] = {'Z': Z, 'list' : formatted_res}
            '''
        
    
    return(kappaClusters)

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




