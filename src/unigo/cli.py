
from .utils import loadUniprotIDsFromCliFiles, unigo_tree_from_api
import json
from .stat_utils import applyOraToVector
from . import uloads as createGOTreeFromAPI

def run(expUniprotIdFile, deltaUniprotIdFile, goApiPort, taxid, method="fisher", asVector=True):
    expUniprotID, deltaUniprotID = loadUniprotIDsFromCliFiles(\
                                            expUniprotIdFile,\
                                            deltaUniprotIdFile
                                            )
    if asVector: # Vector ora     
        print("Running the vectorized ora")
        resp = utils.unigo_vector_from_api(goApiPort, taxid)
        if resp.status_code != 200:
            print(f"request returned {resp.status_code}")
        
        vectorizedProteomeTree = json.loads(resp.text)
        res = applyOraToVector(vectorizedProteomeTree, expUniprotID, deltaUniprotID, 0.5)
        print(res)

    else:# Tree ora       
        resp = utils.unigo_tree_from_api(goApiPort, taxid)
        if resp.status_code != 200:
            print(f"request returned {resp.status_code}")  

        unigoTreeFromAPI = createGOTreeFromAPI(resp.text, expUniprotID)
        x,y = unigoTreeFromAPI.dimensions
        assert not unigoTreeFromAPI.isExpEmpty
        if method == "fisher":
            print("Computing ORA")
            rankingsORA = unigoTreeFromAPI.computeORA(deltaUniprotID)
            print(f"Test Top - {nTop}\n{rankingsORA[:nTop]}")
    