from ...api.store.client import unigoList
from prompt_toolkit import print_formatted_text, HTML
from .connect import bConnect
from . import signatureCheck
from pprint import pformat, pprint
from prompt_toolkit.styles import Style
from ...utils import get_available_uniprot_collection

@bConnect
@signatureCheck
def tlist(*args):
    d = unigoList(*args)
    # The style sheet.
    """
    _style = {
        ':F': '#ff0066',
        ':P': '#ff0066',
        ':C': '#ff0066',
    }
    """
    #style =  Style.from_dict( _style )

    print_formatted_text(HTML("------ <b>Database content</b> --------".center(80)) )
    _ = pformat(d, depth=4, sort_dicts=True, width=80)
    print_formatted_text(_)
    # Print stylinf

def ulist(host, port):
    ucoll_view = get_available_uniprot_collection(host, port)
    print_formatted_text(HTML("------ <b>Available Uniprot collections</b> --------".center(80)) )
    _ = pformat(ucoll_view, depth=4, sort_dicts=True, width=80)
    print_formatted_text(_)