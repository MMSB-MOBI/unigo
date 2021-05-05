from ...api.store.client import unigoList
from prompt_toolkit import print_formatted_text, HTML
from .connect import bConnect
from . import signatureCheck
from pprint import pformat, pprint
from prompt_toolkit.styles import Style


@bConnect
@signatureCheck
def clist(*args):
    d = unigoList(*args)
    # The style sheet.
    _style = {
        ':F': '#ff0066',
        ':P': '#ff0066',
        ':C': '#ff0066',
    }
    
    style =  Style.from_dict( _style )

    print_formatted_text(HTML("------ <b>Database content</b> --------".center(80)) )
    _ = pformat(d, depth=4, sort_dicts=True, width=80)
    print(_)
    # Print stylinf