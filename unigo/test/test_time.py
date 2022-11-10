import sys
sys.path.append('../../')
import unigo

from pydantic import BaseModel
from typing import List

class GODatum(BaseModel):
    id: str
    evidence: str
    term: str

class ProteinObject(BaseModel):
    id: str
    go: List[GODatum]

import json
proteins = json.load(open('../../notebooks/data/proteins.json'))

def get_iterator(json_models):
    for p in json_models:
        yield ProteinObject.parse_obj(p)

protein_iterator = get_iterator(proteins) 

ontology = "/data1/cecile/PSF/ontology/go_2207.owl" 

unigo.setOntology(ontology)
tree = unigo.tree.createGoTree('biological process', protein_iterator)

