from prompt_toolkit import print_formatted_text, HTML
from prompt_toolkit.completion import NestedCompleter
from prompt_toolkit.auto_suggest import AutoSuggest, Suggestion
from decorator import decorator


signatures = {
    'clist': { 
        'vector',
        'tree',
        'culled',
        'all'
    },
    'connect' : None,    
    'load'    : None,
    'exit'    : None,

}
completer = NestedCompleter.from_nested_dict(signatures)

class SignatureBaseError(Exception):
    def __init__(self, cmd):
        self.cmd = cmd
        self.tokens = set(signatures.keys())
        if cmd in signatures:
            self.parameters = signatures[cmd] if type(signatures[cmd]) == set else set(signatures[cmd].keys())
        super().__init__()
        
class SignatureCallError(SignatureBaseError):
    def __init__(self, cmd):
        super().__init__(cmd)#self.message)
    def __str__(self):
        return f"<ansired><u>{self.cmd}</u> is not a valid command call </ansired> <ansigreen><i>{self.tokens}</i><ansigreen>"

class SignatureWrongParamError(SignatureBaseError):
    def __init__(self, cmd, *params):
        self.cmd = cmd
        super().__init__(cmd)#self.message)
        self.badParams = params
    def __str__(self):
        return f"<ansired>{self.cmd} was called with wrong parameters <u>{self.badParams}</u> instead of <i>{self.parameters}</i></ansired>"

class SignatureEmptyParamError(SignatureBaseError):
    def __init__(self, cmd):
        super().__init__(cmd)#self.message)
    def __str__(self):
        return f"<ansired><u>{self.cmd}</u> empty argument list instead of <i>{self.parameters}</i></ansired>"

# Should be made recursive if needed
# We should get function name
@decorator
def signatureCheck(fn, *args, **kwargs):
    cmd = fn.__name__
    if not cmd in signatures:
        raise SignatureCallError(cmd)
    availArg = set([ subCmd for subCmd in signatures[cmd] ])
    _ = set(args) - set(availArg)
    if not args:
        raise SignatureEmptyParamError(cmd)
        #print_formatted_text(HTML(f"<ansired><u>{cmdName}</u> empty argument list</ansired>"))
        #return
    if _:
        raise SignatureWrongParamError(cmd, *_)
        #print_formatted_text( HTML(f"<ansired><u>{cmdName}</u> unkwnow argument <u>{_}</u></ansired>") )
        #return 
    return fn(*args, **kwargs)


class customAutoSuggest(AutoSuggest):
    def __init__(self):
        super().__init__()
    def update(self, cmd):
        self.currentCommand = cmd

    def get_suggestion(self, buffer, document):
        #print("###", document.text)
        if str(document.text).startswith('connect'):
            return Suggestion(" 127.0.0.1 1234")

        #if str(document.text).startswith('list'):
        #    return Suggestion(" all|vector|culled|tree")

        return Suggestion("")