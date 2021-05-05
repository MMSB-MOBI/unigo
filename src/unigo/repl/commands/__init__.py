from decorator import decorator
from prompt_toolkit import print_formatted_text, HTML

from ..helpers import signatureCheck, signatureCheck
from .mutators import load
from .viewers import clist
from .connect import connect, getPrompt

class Executor():
    def __init__(self):
        self.executor = {
            'clist'    : clist,
            'connect' : connect,
            'load'    : load,
            '_exit'    : _exit
        }
    
    def isa(self, cmd):
        return cmd in self.executor
    
    def process(self, cmd, *args):
        if not self.isa(cmd):
            raise ValueError(f"{cmd} is not a valid command {self.executor.keys()}")
        fn = self.executor[cmd]
        
        try:
            fn(*args)
        except Exception as e:
            print_formatted_text(HTML(e))
            return False
        return True

def _exit():
    print_formatted_text(HTML(f"\n\n<skyblue>See you space cowboy</skyblue>"))
    exit(0)

"""
def stripTaxid(data):
    ok     = []
    notok  = []
    for k in data if re.match("^[^\:]")
"""






