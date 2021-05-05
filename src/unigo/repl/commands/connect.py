from decorator import decorator
from ...api.store.client import handshake
from prompt_toolkit import print_formatted_text, HTML

HOST = None
PORT = None
INV_PROMPT = ">>>"

def getPrompt():
    return INV_PROMPT

@decorator
def bConnect(f, *args, **kwargs):
    if HOST is None:
        print_formatted_text(HTML(f"<ansired>Please <u>connect</u> first</ansired>"))
        return
    return f(*args, **kwargs)

def connect(*args):
    global HOST, PORT, INV_PROMPT
    if len(args) != 2:
        print_formatted_text(HTML(f"<ansired><u>connect</u> expects 2 arguments</ansired>"))
        return False
    print(f"Will connect at {args[0]} {args[1]}")
    try : 
        _ = handshake(*args)
        HOST=args[0]
        PORT=args[1]
        print_formatted_text(HTML(f"<ansigreen>Connection successfull at {HOST}:{PORT}</ansigreen>"))
        INV_PROMPT = f"{HOST}:{PORT}>"

        
    except :
        print_formatted_text(HTML(f"<ansired>Connection failed at {args[0]}:{args[1]} </ansired>"))