
from .utils import loadUniprotIDsFromCliFiles, unigo_vector_from_api, unigo_tree_from_api
import json
from .stat_utils import applyOraToVector
from . import uloads as createGOTreeFromAPI

# NO CULLING RESSOURCE USAGE, TO IMPLEMENT
def run(expUniprotIdFile, deltaUniprotIdFile, goApiHost, goApiPort, taxid, method="fisher", asVector=True):
    expUniprotID, deltaUniprotID = loadUniprotIDsFromCliFiles(\
                                            expUniprotIdFile,\
                                            deltaUniprotIdFile
                                            )

    print(expUniprotIdFile, deltaUniprotIdFile, goApiHost, goApiPort, taxid, method)

    if asVector: # Vector ora     
        
        resp = unigo_vector_from_api(goApiHost, goApiPort, taxid)
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

    else:# Tree ora       
        resp = unigo_tree_from_api(goApiHost, goApiPort, taxid)
        if resp.status_code != 200:
            print(f"request returned {resp.status_code}")  

        unigoTreeFromAPI = createGOTreeFromAPI(resp.text, expUniprotID)
        x,y = unigoTreeFromAPI.dimensions
        assert not unigoTreeFromAPI.isExpEmpty
        if method == "fisher":
            print("Computing ORA")
            rankingsORA = unigoTreeFromAPI.computeORA(deltaUniprotID)
            print(f"Test Top - {nTop}\n{rankingsORA[:nTop]}")
    