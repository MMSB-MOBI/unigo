from prompt_toolkit import print_formatted_text, HTML
from prompt_toolkit.completion import NestedCompleter
from prompt_toolkit.auto_suggest import AutoSuggest, Suggestion
from decorator import decorator

import re
"""
 CLI command Signature declaration and check
 we check command and subcommand

 To be able to check parameters, we could declare

 declare = {
     'connect' : { [""]
 }
"""

signatures = {
    'clist': { 
        'vector',
        'tree',
        'culled',
        'all'
    },
    'connect' : None,    
    'build'   : None,
    'load'    : None,
    'exit'    : None,
    'delete'  : None,
    'help'    : {'clist', 'connect',  'load', 'exit', 'delete', 'build'}

}

completer = NestedCompleter.from_nested_dict(signatures)

class SignatureBaseError(Exception):
    def __init__(self, cmd):
        self.cmd = cmd
        self.availbleCmd = set(signatures.keys())
        if cmd in self.availbleCmd:
            self.subCmd = signatures[cmd] if type(signatures[cmd]) == set else set(signatures[cmd].keys())
        super().__init__()
        
class SignatureCallError(SignatureBaseError):
    def __init__(self, cmd):
        super().__init__(cmd)#self.message)
    def __str__(self):
        return f"<ansired><u>{self.cmd}</u> is not a valid command call </ansired> <ansigreen><i>{self.availbleCmd}</i><ansigreen>"

class SignatureWrongSubCommandError(SignatureBaseError):
    def __init__(self, cmd, *subCmd):
        super().__init__(cmd)#self.message)
        self.badSubCmd = subCmd
    def __str__(self):
        return f"<ansired>{self.cmd} was called with wrong subcommand <u>{self.badSubCmd}</u> instead of <i>{self.subCmd}</i></ansired>"

class SignatureEmptySubCommandError(SignatureBaseError):
    def __init__(self, cmd):
        super().__init__(cmd)#self.message)
    def __str__(self):
        return f"<ansired><u>{self.cmd}</u> empty argument list instead of <i>{self.subCmd}</i></ansired>"

# Should be made recursive if needed
@decorator
def signatureCheck(fn, *args, **kwargs):
    cmd = fn.__name__
    if not cmd in signatures:
        raise SignatureCallError(cmd)
    if signatures[cmd] is None:
        return fn(*args, **kwargs)
    availArg = set([ subCmd for subCmd in signatures[cmd] ])
    
    _ = set(args) - set(availArg)
    if not args:
        raise SignatureEmptySubCommandError(cmd)        
    if _:
        raise SignatureWrongSubCommandError(cmd, *_)        
    return fn(*args, **kwargs)

class customAutoSuggest(AutoSuggest):
    def __init__(self):
        super().__init__()
    def update(self, cmd):
        self.currentCommand = cmd

    def get_suggestion(self, buffer, document):
        if str(document.text).startswith('connect'):
            _ = re.findall('([\S]+)', document.text)
            if len(_) == 1:
                return Suggestion(" 127.0.0.1 1234")
            if len(_) == 2:
                return Suggestion(" 1234")
            
        #if str(document.text).startswith('list'):
        #    return Suggestion(" all|vector|culled|tree")
        if str(document.text).startswith('load'):
            _ = re.findall('([\S]+)', document.text)
            if len(_) > 1:
                return Suggestion( pathSuggest(_[-1]) )

        return Suggestion("")

import glob 
from os.path import commonprefix
def pathSuggest(string):
    _ = glob.glob(f"{string}*")
    return re.sub(rf'^{string}', '',commonprefix(_) )
