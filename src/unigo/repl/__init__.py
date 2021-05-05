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
    if not executor.isa(cmd):
        print_formatted_text(HTML(f"<ansired><b>{cmd}</b> is not a valid command</ansired>"))
        return False

    _ = executor.process(cmd, *_)


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
