from marshmallow import Schema, fields, pre_load
from pyproteinsExt.uniprot import isValidID
from ..tree import enumNS
from .. import Univgo as createUnivgo

DEFAULT_MIN_COUNT = 0
DEFAULT_MAX_COUNT = 50
DEFAULT_MAX_FREQ = 1.0

def validateCount(n):
    if n < 0:
        raise ValidationError("Count must be greater than 0.")

def validateFreq(n):
    if n < 0:
        raise ValidationError("maxFreq must be greater than 0.")
    if n > 1.0:
        raise ValidationError("maxFreq must be lower or equal to 1.")
    
class CulledGoParameters:
    def __init__(self, minCount, maxCount, maxFreq):
        self.minCount = int(minCount)
        self.maxCount = int(maxCount)
        self.maxFreq  = float(maxFreq)
    def __repr__(self):
        return f'<CulledGoParameters(minCount={self.minCount}, maxCount={self.maxCount}, maxFreq={self.maxFreq})>' 

class CulledGoParametersSchema(Schema):
    global DEFAULT_MIN_COUNT, DEFAULT_MAX_COUNT, DEFAULT_MAX_FREQ

    @pre_load(pass_many=True)
    def unwrap_request(self, request, **kwargs):
        return request.get_json()

    minCount = fields.Int(validate=validateCount,  missing=DEFAULT_MIN_COUNT)
    maxCount = fields.Int(validate=validateCount,  missing=DEFAULT_MAX_COUNT)
    maxFreq  = fields.Float(validate=validateFreq, missing=DEFAULT_MAX_FREQ)

def validateUnivgoSerial(fn):
    
    def loadUnivGO_dec(data):
        for k in ["single_tree", "omega_uniprotID", "ns"]:
            if k not in data:
                raise KeyError(f"univoGo serial key {k} missing")
    
        if not type(data["omega_uniprotID"]) == list:
            raise KeyError(f"univGo serial key omega_uniprotID is not a list\n=>{type(data['omega_uniprotID'])}")
        for uID in data["omega_uniprotID"]:
            if not isValidID(uID):
                raise ValueError(f"univGo serial uniprotID not valid {uID}")
        if data["ns"] not in enumNS:
            raise ValueError(f"univGo serial NS not valid {data['ns']}")
        return fn(data)
        
    
    return loadUnivGO_dec

@validateUnivgoSerial
def loadUnivGO(d):
    univGo = createUnivgo(serial=d)
    return univGo
