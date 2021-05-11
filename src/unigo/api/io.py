from flask import request, abort
from .. import utils
def checkPwasInput(enforcedCulling=True):
    data = request.get_json()

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