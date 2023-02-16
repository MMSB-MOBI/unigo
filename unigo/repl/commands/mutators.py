from prompt_toolkit import print_formatted_text, HTML
from prompt_toolkit.shortcuts import ProgressBar
from .connect import bConnect
from . import signatureCheck
from ...utils.uniprot import sync_to_uniprot_store
from ...api.store.client import addTree3NSByTaxid as goStoreAdd, delTaxonomy, buildVectors
from ...api.store.client import InsertionError, DeletionError
import time

@bConnect
@signatureCheck
def load(owlFile, *args):
    print(f"Connecting to goStore service to add following proteome(s) {args}")
    
    try:
        taxidTreeIter = sync_to_uniprot_store("localhost", 6379, owlFile, *args)
        msg = goStoreAdd(taxidTreeIter, fromCli=True)
        print_formatted_text(HTML(f'<ansigreen>{msg}</ansigreen>'))
    except InsertionError as e:
        print_formatted_text(HTML(f'<ansired>{e}</ansired>'))
    except KeyError as e:
        print_formatted_text(HTML(f'<ansired>{e}</ansired>'))

@bConnect
@signatureCheck
def delete(*taxids):
    try:
        msg = delTaxonomy(taxids, fromCli=True)
    except DeletionError as e:
        print_formatted_text(HTML(f'<ansired>{e}</ansired>'))
    
    print_formatted_text(HTML(f'<ansigreen>{msg}</ansigreen>'))

@bConnect
@signatureCheck
def build():
    try:
        bVectorStatusChange = buildVectors(fromCli=True)
        status, size = next(bVectorStatusChange)
        if status == "nothing to build":
            print_formatted_text( HTML(f'<ansigreen>No vector to build</ansigreen>') )
            return

        if status == "running":
            print_formatted_text( HTML(f'<ansigreen>A building thread is already running, monitoring it...</ansigreen>') )
        else:
            print_formatted_text( HTML(f'<ansigreen>Start Building</ansigreen>')) 
        with ProgressBar() as pb:
            for i in pb( bVectorStatusChange, total=int(size) ):
                time.sleep(0.1)
    except Exception  as e:
        print_formatted_text( HTML(f'<ansired>Error while Building {e}</ansired>') )