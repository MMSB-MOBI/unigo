from prompt_toolkit import PromptSession, print_formatted_text, HTML
from prompt_toolkit.key_binding import KeyBindings
from .helpers import completer, customAutoSuggest
import re
from .commands import getPrompt, Executor


executor = Executor()

def digest(_input):
    if re.match('^[\s]*$', _input) :
        #print("Empty string")
        return
    _ = _input.split()
    #print("input::", _)
    
    cmd = _.pop(0)
    if cmd == 'help':
        _help(_[0] if _ else None)
        return True
    if not executor.isa(cmd):
        print_formatted_text(HTML(f"<ansired><b>{cmd}</b> is not a valid command</ansired>\nType help <command> for details"))
        return False

    _ = executor.process(cmd, *_)

def _help(cmd):
    if cmd is None:
        msg = ', '.join([ f"<u>{c}</u>" for c in executor.availbleCmd ])
        msg = f"Available command: {msg}"
        print_formatted_text(HTML(msg))
    elif not executor.isa(cmd):
        print_formatted_text(HTML(f"<ansired><u>{cmd}</u> is not a valid command</ansired>\n{msg}"))
    else:
        msg = executor.help(cmd)
        print_formatted_text(HTML(f"<u>{cmd}</u>:\n{msg}"))


kb = KeyBindings()
@kb.add('c-c')
def exit_(event):
    """
    Pressing Ctrl-Q will exit the user interface.

    """
    #_exit()
    executor.process('exit')


def run():
    session = PromptSession()
    
    while True:
        answer = session.prompt(getPrompt(), completer=completer, auto_suggest=customAutoSuggest(),
        key_bindings=kb )
        digest(answer)
        #print('You said: %s' % answer)
    exit(1)
