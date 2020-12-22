"""Go ontology tree manipulation tool and microservice

Usage:
    unigo run (tree|fisher|convert) [--onto=<owlFile>] [--prot=<xmlFile>] [--size=<n_exp_proteins>] [--silent] [--delta=<n_modified_proteins>] [--head=<n_best_pvalue>]
    unigo api [--port=<portNumber>] [--onto=<owlFile>] [--prot=<xmlFile>]
    unigo client (tree|fisher) <taxid> [--port=<portNumber>] [--silent]

Options:
  -h --help     Show this screen.
  --onto=<owlFile> ontology file location in owl format [default:<package_location>/data/go.owl]
  --prot=<xmlFile> uniprot file location in xml format [default:<package_location>/data/uniprot-proteome_UP000000807.xml.gz]
  --size=<number> number of uniprotID to build the experimental protein group from the uniprot elements [default:50]
  --delta=<number> fraction of the experimental protein group to build the group of over/undee-represented protein [default:0.1]
  --silent  stop ORA scoring dump
  --head=<n_best_pvalue> display n best GO pathway [default:5]

"""
import os

DEFAULT_PROTEOME=f"{os.path.dirname(os.path.abspath(__file__))}/data/uniprot-proteome_UP000000807.xml.gz"

from docopt import docopt
from pyproteinsExt import uniprot as pExt
from . import Unigo as createGOTreeTest
from . import uloads as createGOTreeTestFromAPI

from . import Univgo as createGOTreeTestUniverse
from . import vloads as createGOTreeTestUniverseFromAPI

from .api import listen
from requests import get

arguments = docopt(__doc__)
print(arguments)

nDummy = int(arguments['--size']) if arguments['--size'] else 50
nTop   = int(arguments['--head']) if arguments['--head'] else 5
proteomeXML = arguments['--prot'] if arguments['--prot'] else DEFAULT_PROTEOME
apiPort = arguments['--port'] if arguments['--port'] else 5000



if arguments['fisher'] or arguments['convert']:
    arguments['tree'] = True


if arguments['run'] or arguments['api'] or arguments['client']:
    uColl = pExt.EntrySet(collectionXML=proteomeXML)
    print(f"Setting up a dummy experimental collection of {nDummy} elements")
    expUniprotID =[]
    for uObj in uColl:
        if len(expUniprotID) == nDummy:
            break
        if uObj.isGOannot:
            expUniprotID.append(uObj.id)
    
    nDelta = int(nDummy*0.1)
    print(f"Considering {nDelta} proteins among {nDummy} experimental as of significantly modified quantities")
    deltaUniprotID = expUniprotID[:nDelta]

if arguments['api']:
    print("Testing universe tree API")
    print(f"Loading Test proteome{proteomeXML} as universe")
    uColl = pExt.EntrySet(collectionXML=proteomeXML)
    tree_universe = createGOTreeTestUniverse( 
                                    ns          = "biological process", 
                                    fetchLatest = False,
                                    uniColl     = uColl)

    app = listen( trees=[tree_universe], taxids=[uColl.taxids] )
    app.run(debug=False)


if arguments['run']:
    print("Testing local implementation")
    if arguments['tree']:
        print("Creatin unigo Tree")
        
        unigoTree = createGOTreeTest(backgroundUniColl = uColl,
                                    proteinList       = expUniprotID,
                                    fetchLatest       = False)
        
    if arguments['fisher']:
        print("Computing ORA")
        rankingsORA = unigoTree.computeORA(deltaUniprotID, verbose = not arguments['--silent'])
        print(f"Test Top - {nTop}\n{rankingsORA[:nTop]}")

    if arguments['convert']:
        print("Testing tree serialization")
        d = unigoTree.tree.f_serialize()
        print(d.asDict)



if arguments['client']:
    url = f"http://127.0.0.1:{apiPort}/unigo/{arguments['<taxid>']}"
    print(f"Testing Annotation tree API from {url}")
    
    resp = get(url)
    if resp.status_code == 404:
        print(f"{url} returned 404, your taxid {arguments['<taxid>']} may not be registred")
    else:
        print(f"Trying to build a single universal based on taxid {arguments['<taxid>']} API response")
        tree_universe = createGOTreeTestUniverseFromAPI(resp.text)
        n, l, p, p_nr = tree_universe.dimensions
        print(f"Successfully loaded a single {n} nodes and {l} links universal tree")
        print("Trying to create Unigo Object using previously set experimental protein list")
        unigoTreeFromAPI = createGOTreeTestFromAPI(resp.text, expUniprotID)
        x,y = unigoTreeFromAPI.dimensions
        print("Unigo Object successfully buildt w/ following dimensions:")
        print(f"\txpTree => nodes:{x[0]} children_links:{x[1]}, total_protein_occurences:{x[2]}, protein_set:{x[3]}")  
        print(f"\t universeTree => nodes:{y[0]} children_links:{y[1]}, total_protein_occurences:{y[2]}, protein_set:{y[3]}")  
        if arguments['fisher']:
            print("Computing ORA")
            rankingsORA = unigoTreeFromAPI.computeORA(deltaUniprotID, verbose = not arguments['--silent'])
            print(f"Test Top - {nTop}\n{rankingsORA[:nTop]}")