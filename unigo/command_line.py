
from .api.store.client.viewers import unigo_vector_from_api, unigo_tree_from_api
from .utils.io import loadUniprotIDsFromCliFiles
import json
from .stats.ora import applyOraToVector
from . import Unigo

"""
WIP local run
"""
def run(expUniprotIdFile, deltaUniprotIdFile, taxid, method="fisher", asVector=True, n_top=10):
    expUniprotID, deltaUniprotID = loadUniprotIDsFromCliFiles(\
                                            expUniprotIdFile,\
                                            deltaUniprotIdFile
                                            )
  
    if asVector: # Vector ora     
        
        resp = unigo_vector_from_api(taxid)
        if resp.status_code != 200:
            print(f"request returned {resp.status_code}")
            return None
        ns='molecular function'
        print(f"Running the vectorized ora on {ns}")
        vectorizedProteomeTree = json.loads(resp.text)
       # print(vectorizedProteomeTree[ns].keys())
       # print( len(expUniprotID), len(deltaUniprotID) )
        res = applyOraToVector(vectorizedProteomeTree[ns], expUniprotID, deltaUniprotID, 0.05, translateID=True)
        print(res)

    else:# Tree ora       TO DO AFTER TEST INSTANCE
        resp = unigo_tree_from_api(taxid)
        if resp.status_code != 200:
            print(f"request returned {resp.status_code}")  
        # desrialize blueprint
        unigo_blueprint = Unigo(from_serial = resp.json()) 
        #unigoTreeFromAPI = createGOTreeFromAPI(resp.text, expUniprotID)
        x,y = unigo_blueprint.dimensions

        #assert not unigo_blueprint.isExpEmpty
        if method == "fisher":
            print("Computing ORA")
            # uncomment and fix below
            #rankingsORA = unigoTreeFromAPI.computeORA(deltaUniprotID)
            #print(f"Test Top - {n_top}\n{rankingsORA[:n_top]}")
    