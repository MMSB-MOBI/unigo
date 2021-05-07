from decorator import decorator
from prompt_toolkit import print_formatted_text, HTML

from ..helpers import signatureCheck, signatureCheck
from .mutators import load, delete, build
from .viewers import clist
from .connect import connect, getPrompt
from decorator import decorator
"""
class ExecutorBaseError(Exception):
    def __init__(self, symbol):
        self.typings = cmd
        self.availbleCmd = set(signatures.keys())
        if cmd in availbleCmd:
            self.subCmd = signatures[cmd] if type(signatures[cmd]) == set else set(signatures[cmd].keys())
        super().__init__()
        
class ExecutorParameterNumberError(SignatureBaseError):
    def __init__(self, cmd):
        super().__init__(cmd)#self.message)
    def __str__(self):
        return f"<ansired><u>{self.cmd}</u> is not a valid command call </ansired> <ansigreen><i>{self.availbleCmd}</i><ansigreen>"

class ExecutorParameterNumberError(SignatureBaseError):
    def __init__(self, cmd):
        super().__init__(cmd)#self.message)
    def __str__(self):
        return f"<ansired><u>{self.cmd}</u> is not a valid command call </ansired> <ansigreen><i>{self.availbleCmd}</i><ansigreen>"
"""



class Executor():
    def __init__(self):
        self.executor = {
            "clist" : {
                "_target" : clist,
                "paramTypes" : None,
                "help" : f"List database content\nUsage:clist tree|all|culled|vector"
            },
            "connect" : {
                "_target"    : connect,
                "paramTypes" : [str, int],
                "runningParams" : True,
                "help" : f"Connect to database\nUsage: connect <i>hostname</i> <i>port</i>"
            },
            "load"    :  {
                "_target"    : load,
                "paramTypes" : [str, str],
                "runningParams" : True,
                "help" : f"Insert any number of triplets of GO trees in the database\n Usage:load <i>owl_file proteome_file_1, ...</i>"
            },
            'exit'    :  {
                "_target"    : _exit,
                "paramTypes" : None,
                "help" : f"Close the interface\nUsage: exit"
            },
            "delete"    :  {
                "_target"    : delete,
                "paramTypes" : [str],
                "runningParams" : True,
                "help" : f"Delete any number of triplets of GO trees in the database\nUsage: delete <i>taxid_1, ...</i>"
            },
            "build"    :  {
                "_target"    : build,
                "paramTypes" : None,               
                "help" : f"Build all missing triplets of <u>GO vectors</u> in the database\nUsage: build</i>"
            }
        }
    
    @property
    def availbleCmd(self):
        for k in self.executor:
            yield k
        
    def isa(self, cmd):
        return cmd in self.executor
    
    def help(self, cmd):
        return self.executor[cmd]["help"]

    def process(self, cmd, *args):
        if not self.isa(cmd):
            raise ValueError(f"{cmd} is not a valid command {self.executor.keys()}")
        fn = self.executor[cmd]["_target"]
        try:
            self.checkTypes(cmd, *args)
            fn(*args)
        except Exception as e:
            print_formatted_text(HTML(e))
            return False
        return True

    def checkTypes(self, symbol, *args):
        _ = self.executor[symbol]
        typeArgs = _["paramTypes"]
        if typeArgs is None:
            return True
        
        if len(args) < len(typeArgs) or\
            (len(args) > len(typeArgs) and _["runningParams"]) :
            raise TypeError(f"<ansired>{symbol} excepts {len(typeArgs)} argument(s)</ansired>")
""" GIVING UP ON PARAM TYPE CHECKING...
        for i, carg in enumerate(args):
            ctype=type(carg)
            j = i if i < len(typeArgs) else len(typeArgs) - 1
            _type = typeArgs[j]
            
            if ctype != _type:
                _type = str(_type).replace('<', '').replace('>', '').replace("'", "")
                raise TypeError(f"<ansired>{symbol} argument {carg} wrongs types should be {_type}</ansired>")
                #print(f"<ansired>{symbol} argument {carg} has wrong {ctype}</ansired>")
                #raise TypeError(f"<ansired>{symbol} argument {carg} has wrong {ctype}</ansired>")
"""
def _exit():
    print_formatted_text(HTML(f"\n\n<skyblue>See you space cowboy</skyblue>"))
    exit(0)

"""
def stripTaxid(data):
    ok     = []
    notok  = []
    for k in data if re.match("^[^\:]")
"""






