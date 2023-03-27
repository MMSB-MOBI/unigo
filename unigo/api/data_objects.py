from marshmallow import Schema, fields, pre_load
from pyproteinsext.uniprot import isValidID
from ..tree import enumNS
from .. import Unigo


DEFAULT_MIN_COUNT = 0
DEFAULT_MAX_COUNT = 50
DEFAULT_MAX_FREQ  = 1.0
DEFAULT_PVALUE    = 0.05

class ValidationError(Exception):
    pass

def validateCount(n):
    if n < 0:
        raise ValidationError("Count must be greater than 0.")

def validate_pvalue(n):
    if n < 0 or n > 1:
        raise ValidationError("pvalue must be between than 0. and 1.")

def validateFreq(n):
    if n < 0:
        raise ValidationError("maxFreq must be greater than 0.")
    if n > 1.0:
        raise ValidationError("maxFreq must be lower or equal to 1.")
    
class CulledGoParameters:
    def __init__(self, minCount, maxCount, maxFreq, pvalue):
        self.minCount = int(minCount)
        self.maxCount = int(maxCount)
        self.maxFreq  = float(maxFreq)
        self.pvalue  = float(pvalue)
    def __repr__(self):
        return f'<CulledGoParameters(minCount={self.minCount}, maxCount={self.maxCount}, maxFreq={self.maxFreq}, pvalue={self.pvalue})>' 

class CulledGoParametersSchema(Schema):
    global DEFAULT_MIN_COUNT, DEFAULT_MAX_COUNT, DEFAULT_MAX_FREQ, DEFAULT_PVALUE

    @pre_load(pass_many=True)
    def unwrap_request(self, request, **kwargs):
        return request.get_json()

    minCount = fields.Int(validate=validateCount,  missing=DEFAULT_MIN_COUNT)
    maxCount = fields.Int(validate=validateCount,  missing=DEFAULT_MAX_COUNT)
    maxFreq  = fields.Float(validate=validateFreq, missing=DEFAULT_MAX_FREQ)
    pvalue   = fields.Float(validate=validate_pvalue, missing=DEFAULT_PVALUE)

def validate_unigo_serial(fn):
    
    def load_unigo_dec(data):
        for k in ["tree", "omega_uniprotID", "ns"]:
            if k not in data:
                raise KeyError(f"unigo serial key {k} missing from [{list(data.keys())}]")
    
        if not type(data["omega_uniprotID"]) == list:
            raise KeyError(f"unigo serial key omega_uniprotID is not a list\n=>{type(data['omega_uniprotID'])}")
        for uID in data["omega_uniprotID"]:
            if not isValidID(uID):
                raise ValueError(f"unigo serial uniprotID not valid {uID}")
        if data["ns"] not in enumNS:
            raise ValueError(f"unigo serial NS not valid {data['ns']}")
        return fn(data)
        
    
    return load_unigo_dec

@validate_unigo_serial
def unigo_deserializer(d):
    my_unigo = Unigo(from_serial=d)
    return my_unigo
