
from unigo.api.io import check_pwas_input_from_cmdline
from unigo.api.pwas_compute import compute_over_vector
from .utils import loadUniprotIDsFromCliFiles, unigo_vector_from_api, unigo_tree_from_api
import json
from .stat_utils import applyOraToVector
from . import uloads as createGOTreeFromAPI
from .api.io import check_pwas_input_from_cmdline as check_pwas_input

# NO CULLING RESSOURCE USAGE, TO IMPLEMENT
def _run(expUniprotIdFile, deltaUniprotIdFile, goApiHost, goApiPort, taxid, method="fisher", asVector=True):
    expUniprotID, deltaUniprotID = loadUniprotIDsFromCliFiles(\
                                            expUniprotIdFile,\
                                            deltaUniprotIdFile
                                            )

    #print(expUniprotIdFile, deltaUniprotIdFile, goApiHost, goApiPort, taxid, method)

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
        #print(res)

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

def run(exp_uniprot_id_file, delta_uniprot_id_file, go_host, go_port, taxid, min_count = 0, max_count = 50, max_freq = 1, method="fisher", compute_clusters = True, pvalue=0.05, force_universal=False):
    exp_ids, delta_ids = loadUniprotIDsFromCliFiles(\
                                            exp_uniprot_id_file,\
                                            delta_uniprot_id_file
                                            )
    data = check_pwas_input(exp_ids, delta_ids, taxid, method)
    data['minCount'] = min_count
    data['maxCount'] = max_count
    data['maxFreq'] = max_freq
    data['pvalue'] = pvalue
    results = compute_over_vector(data, go_host, go_port, compute_clusters, force_universal=force_universal)
    return results
    
