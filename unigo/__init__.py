from .tree import setOntology, createGoTree, load, enumNS
import os
import json

import datetime

print("DVL PYPROTEINS_EXT PKG", datetime.datetime.now())

"""
Minimal ontology tree of a uniprot collection

Merging Uniprot Object Collection 
Consrtuctor

"""

class Unigo:
    def __init__(self, from_prev=None,
                       from_serial = None, 
                       from_json_string = None,
                       uniColl = None, 
                       owlFile=None, fetchLatest=True,
                       ns="biological process"):#, **kwargs):
        self.masked = False

        self._GO_index = {}

        if not from_prev is None:           
            ns, tree = from_prev
            self.tree = tree
            self.ns = ns
            self.omega_uniprotID = list(set(self.tree.proteins))
            return

        if not ns in enumNS:
            raise TypeError(f"Provided namespace {ns} is not registred in {enumNS}")

        if from_json_string or from_serial:
            serial = from_serial if from_serial else json.loads(from_json_string)
            self.tree     = load(serial["tree"])
            self.omega_uniprotID = serial["omega_uniprotID"]
            self.ns              = serial["ns"]
            if not self.ns in enumNS:
                raise TypeError(f"Provided namespace {self.ns} is not registred in {enumNS}")
            self.tree._index()
            return
        
        if uniColl is None :
            raise ValueError(f"parameters uniColl required")
        try :
            if owlFile is None and fetchLatest:
                print("Fetching ontology")
                setOntology(url="http://purl.obolibrary.org/obo/go.owl")#, **kwargs)
            elif not owlFile is None:
                setOntology(owlFile=owlFile)
            else:
                print("Using package release ontology")
                setOntology(owlFile=f"{os.path.dirname(os.path.abspath(__file__))}/data/go.owl")
        except Exception as e:
            print(f"Could not create ontology")
            print(e)
            raise TypeError("Failed creating Go tree")
        # De novo construction, meant to be a blueprint instance, we hence compute background frequencies
        self.tree = createGoTree(        ns = ns,                                  
                                  protein_iterator = uniColl)
        self.tree.compute_background_frequency()
        self.ns = ns
        self.omega_uniprotID = list(set(self.tree.proteins))
        

    def serialize(self):
        return {
            "tree"     : self.tree.f_serialize().asDict,
            "omega_uniprotID" : self.omega_uniprotID,
            "ns"              : self.ns    
        }

    #def picklify(self):
    #    return {
    #        "tree"            : self.tree.makePickable(),
    #        "omega_uniprotID" : self.omega_uniprotID,
    #        "ns"              : self.ns    
    #}



    def vectorize(self):
        data = self.tree.vectorize()
        data["annotated"] = data["registry"][:]
        #data["nb_annot"] = len(data["registry"])
        for uniprotID in set(self.omega_uniprotID) - set(data["registry"]):
            data["registry"].append(uniprotID)
        return data

    @property
    def dimensions(self):
        return self.tree.dimensions
    
    def walk(self, type=None):
        """Iterate over the tree, by default the experimental one"""
        if type is None:
            for node in self.tree.walk():
                yield node

    def getByID(self, id):
        if not self._GO_index:
            self._index_by_GO()
        
        return self._GO_index[id]

    def _index_by_GO(self):
        print("Create GO index")
        for n in self.tree.walk():
            self._GO_index[n.ID] = n
    
"""
Register the location of the provided protein list accross the blueprint tree "masking it"
The tree can then be used to compute ora
Each node in blueprint state si updated w/:
    - etag atttribute's elements now refering to observed uniprot ids
    - background_members attribute now stores the MEMBERS of the annotation node in the blueprint state
Based on the updated eTag elements
"""
def unigo_obs_mask(unigo_blueprint, obs_uniprot_ids)->Unigo:
    if unigo_blueprint.masked:
        raise TypeError("Masking an already masked unigo tree is not allowed")
    
    # Set predicate function to retain node in the blueprint tree which are part of provided proteins
    def keep_proteinID(node):
        return set(node.getMembers()) & set(obs_uniprot_ids)
    # Obtain the tree
    mask_obs_tree = unigo_blueprint.tree.drop(keep_proteinID, noCollapse=True, noLeafCountUpdate=True)
    # Update eTags to keep only provided proteinIDList, and memo member proteins in bkg pop
    for n in mask_obs_tree.walk():
        n_bp = unigo_blueprint.tree.getByID(n.ID) # Get the node in the blueprint topology, as nodes may have been droped in the mask_obs
        n.background_members = n_bp.getMembers(nr=True)
        n.eTag = list(  set(obs_uniprot_ids) & set(n.eTag) )
    # update tree topology and leaf counts
    mask_obs_tree = mask_obs_tree.collapse()
    mask_obs_tree.leafCountUpdate()
    masked_unigo = Unigo(from_prev = (unigo_blueprint.ns, mask_obs_tree))
    masked_unigo.masked = True
    return masked_unigo

def unigo_prune(unigo_obj, predicate):
    """Create a Unigo object by Pruning the supplied Unigo instance

        Parameters
        ----------
        unigoObj : 
        predicate: predicate function applied to prodived Unigo to drop or keep its elements
    """
    tree          = unigo_obj.tree.drop(predicate)
    

    return Unigo( from_prev=(unigo_obj.ns, tree) )


def chop(unigo_obj, name=None, ID=None):
    """Create a Unigo object with its root extracted from the provided Unigo instance

        Parameters
        ----------
        unigoObj : 
        name: name of the node the look for in provided unigoObj. It will be the root of the returned Unigo
        ID: ID of the node the look for in provided unigoObj. It will be the root of the returned Unigo
    """
    
    tree          = unigo_obj.tree.newRoot(name=name, ID=name)
    return Unigo( from_prev=(unigo_obj.ns, tree) )