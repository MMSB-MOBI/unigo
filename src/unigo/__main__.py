"""Go ontology tree manipulation tool and microservice

Usage:
    unigo local_test (tree|fisher|convert) [--onto=<owlFile>] [--prot=<xmlFile>] [--size=<n_exp_proteins>] [--silent] [--delta=<n_modified_proteins>] [--head=<n_best_pvalue>]
    unigo gostore [--goport=<portNumber>] [--onto=<owlFile>] [--prot=<xmlFile>]
    unigo pwas test [--port=<portNumber>] [--method=<statMethod>] [--prot=<xmlFile>] [--silent]
    unigo pwas api [--pwasport=<portNumber>] [--goport=<portNumber>]
    unigo pwas cli --taxid=<number> --exp-prot=<list.txt> --de-prot=<list.txt> [--goport=<portNumber> --method=<statMethod> --vectorized]

Options:
  -h --help     Show this screen.
  --onto=<owlFile> ontology file location in owl format [default:<package_location>/data/go.owl]
  --prot=<xmlFile> uniprot file location in xml format [default:<package_location>/data/uniprot-proteome_UP000000807.xml.gz]
  --size=<number> number of uniprotID to build the experimental protein group from the uniprot elements [default:50]
  --delta=<number> fraction of the experimental protein group to build the group of over/undee-represented protein [default:0.1]
  --silent  stop ORA scoring dump
  --head=<n_best_pvalue> display n best GO pathway [default:5]
  --goport=<portNumber> : port for GO API
  --pwasport=<portNumber> : port for pwas API
  --taxid=<number> : proteome taxid
  --method=<statMethod> : statistical method to compute pathway p-value
  --exp-prot=<list.txt> : txt file with all proteomics experience proteins accession (one per line)
  --de-prot=<list.txt> : txt file with all differentially expressed proteins accession (one per line)
  --vectorized: use the goStore vectorize protocol
"""
import os

DEFAULT_PROTEOME=f"{os.path.dirname(os.path.abspath(__file__))}/data/uniprot-proteome_UP000000807.xml.gz"
DEFAULT_TAXID=243273
DEFAULT_METHOD = "fisher"

from docopt import docopt
from pyproteinsExt import uniprot as pExt
from . import Unigo as createGOTreeTest
from . import uloads as createGOTreeTestFromAPI

from . import Univgo as createGOTreeTestUniverse
from . import vloads as createGOTreeTestUniverseFromAPI

from .api import listen as go_ress_listen
from requests import get

from .stat_utils import applyOraToVector
from . import utils

from .pwas import listen as pwas_listen

arguments = docopt(__doc__)
print(arguments)

nDummy = int(arguments['--size']) if arguments['--size'] else 50
nTop   = int(arguments['--head']) if arguments['--head'] else 5
proteomeXML = arguments['--prot'] if arguments['--prot'] else DEFAULT_PROTEOME
goApiPort = arguments['--goport'] if arguments['--goport'] else 5000
pwasApiPort = arguments['--pwasport'] if arguments["--pwasport"] else 5001
owlFile = arguments['--onto']     if arguments['--onto'] else None
method = arguments['--method'] if arguments['--method'] else DEFAULT_METHOD
taxid = arguments['--taxid'] if arguments["--taxid"] else DEFAULT_TAXID

if arguments['fisher'] or arguments['convert']:
    arguments['tree'] = True

#Load or create proteins sets
if arguments['test']:
    uColl = pExt.EntrySet(collectionXML=proteomeXML)
    print(f"Setting up a dummy experimental collection of {nDummy} elements from {proteomeXML}")
    expUniprotID =[]
    _ = uColl.taxids
    if len(_) > 1:
        print(f"Warnings found many taxids in uniprot collection : {_}")
    taxid = _[0]
    for uObj in uColl:
        if len(expUniprotID) == nDummy:
            break
        if uObj.isGOannot:
            expUniprotID.append(uObj.id)
    
    nDelta = int(nDummy*0.1)
    print(f"Considering {nDelta} proteins among {nDummy} experimental as of significantly modified quantities")
    deltaUniprotID = expUniprotID[:nDelta]

if arguments['cli']:
    expUniprotID = []
    deltaUniprotID = []
    with open(arguments["--exp-prot"]) as f:
        for l in f:
            expUniprotID.append(l.rstrip())
    with open(arguments["--de-prot"]) as f:
        for l in f:
            deltaUniprotID.append(l.rstrip())

    if not utils.check_proteins_subset(expUniprotID, deltaUniprotID):
        raise Exception("Differentially expressed proteins are not completely included in total proteins")

if arguments['gostore']:
    print("Testing universe tree API")
    print(f"Loading Test proteome{proteomeXML} as universe")
    
    try :
        uColl = pExt.EntrySet(collectionXML=proteomeXML)
        tree_universe = createGOTreeTestUniverse( 
                                        owlFile     = owlFile,
                                        ns          = "biological process", 
                                        fetchLatest = False,
                                        uniColl     = uColl)
    except:
        print("Fatal error, early Exit")
        exit(1)
    app = go_ress_listen( trees=[tree_universe], taxids=[uColl.taxids] )
    app.run(debug=False, port=goApiPort)

if arguments['local_test']:
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

if arguments['pwas']:
    if arguments['api']:
        pwas_app = pwas_listen(goApiPort)
        pwas_app.run(debug=True, port=pwasApiPort)
        
    elif not arguments['--vectorized']:         
        resp = utils.unigo_tree_from_api(goApiPort, taxid)
        if resp.status_code != 200:
            print(f"request returned {resp.status_code}")  

        unigoTreeFromAPI = createGOTreeTestFromAPI(resp.text, expUniprotID)
        x,y = unigoTreeFromAPI.dimensions
        assert not unigoTreeFromAPI.isExpEmty
        if method == "fisher":
            print("Computing ORA")
            rankingsORA = unigoTreeFromAPI.computeORA(deltaUniprotID, verbose = not arguments['--silent'])
            print(f"Test Top - {nTop}\n{rankingsORA[:nTop]}")
    else:
        import json
        print("Running the vectorized ora")
        resp = utils.unigo_vector_from_api(goApiPort, taxid)
        if resp.status_code != 200:
            print(f"request returned {resp.status_code}")
        
        vectorizedProteomeTree = json.loads(resp.text)
        res = applyOraToVector(vectorizedProteomeTree, expUniprotID, deltaUniprotID, 0.5)
        print(res)
