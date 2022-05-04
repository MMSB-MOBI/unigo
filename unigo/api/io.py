from typing import TypedDict
from flask import request, abort
from .. import utils

class PwasData(TypedDict):
    all_accessions: list[str]
    taxid: int
    significative_accessions: list[str]
    method: str

def check_pwas_input_from_route():
    data = request.get_json()
    return _check_pwas_input(data)

def check_pwas_input_from_cmdline(exp_uniprot_ids: list[str], delta_uniprot_ids: list[str], taxid: int, method: str):
    data  = {
        'all_accessions' : exp_uniprot_ids,
        'significative_accessions' : delta_uniprot_ids,
        'taxid' : taxid,
        'method' : method
    }

    return _check_pwas_input(data)

def _check_pwas_input(data : PwasData, enforcedCulling=True):
    if not data:
        print(f"ERROR in request : empty posted data. Abort 400")
        abort(400)
    
    for key in ["all_accessions", "taxid", "significative_accessions", "method"]:
        if not key in data:
            print(f"ERROR in request : {key} not in posted data. Abort 400")
            abort(400)

    if not data["method"] in ["fisher"]:
        print(f"ERROR : this statistical method ({data['method']}) is not handled. Availables : fisher")
        abort(400)        

    if not utils.check_proteins_subset(data["all_accessions"], data["significative_accessions"]):
        print(f"ERROR : significative accessions are not all included in all accessions")
        abort(400)
        
    return data



