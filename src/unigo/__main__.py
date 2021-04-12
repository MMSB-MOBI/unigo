"""Go ontology tree manipulation tool and microservice

Usage:
    unigo local_test (tree|fisher|convert) [--onto=<owlFile>] [--prot=<xmlFile>] [--size=<n_exp_proteins>] [--silent] [--delta=<n_modified_proteins>] [--head=<n_best_pvalue>]
    unigo store (start|wipe) [ dry | [ --onto=<owlFile> --prot=<xmlFile> ] ] [--goport=<portNumber>] [--redis --rp=<redis_port> --rh=<redis_host>]
    unigo store (add|del) [--goport=<portNumber> --gohost=<hostname>] [--onto=<owlFile>] [--prot=<xmlFile>]
    unigo pwas test [--port=<portNumber>] [--method=<statMethod>] [--prot=<xmlFile>] [--silent]
    unigo pwas api [--pwasport=<portNumber> --goport=<portNumber>] [ --vectorized]
    unigo pwas cli --taxid=<number> --exp-prot=<list.txt> --de-prot=<list.txt> [--goport=<portNumber> --method=<statMethod> --vectorized]

Options:
  -h --help     Show this screen.
  --onto=<owlFile>  ontology file location in owl format [default:<package_location>/data/go.owl]
  --prot=<xmlFile>  uniprot file location in xml format [default:<package_location>/data/uniprot-proteome_UP000000807.xml.gz]
  --size=<number>  number of uniprotID to build the experimental protein group from the uniprot elements [default: 50]
  --delta=<number>  fraction of the experimental protein group to build the group of over/undee-represented protein [default: 0.1]
  --silent  stop ORA scoring dump
  --head=<n_best_pvalue>  display n best GO pathway [default: 5]
  --goport=<portNumber>  port for GO API [default: 1234]
  --gohost=<hostname>  host name for GO API [default: localhost]
  --pwasport=<portNumber>  port for pwas API [default: 5000].
  --taxid=<number>  proteome taxid
  --method=<statMethod>  statistical method to compute pathway p-value [default: fisher]
  --exp-prot=<list.txt>  txt file with all proteomics experience proteins accession (one per line)
  --de-prot=<list.txt>  txt file with all differentially expressed proteins accession (one per line)
  --vectorized  use the goStore vectorize protocol
  --redis  use redis backend as storage engine
  --rp=<redis_port>  redis DB TCP port [default: 6379]
  --rh=<redis_host>  redis DB http adress [default: localhost]
"""
import os

DEFAULT_PROTEOME=f"{os.path.dirname(os.path.abspath(__file__))}/data/uniprot-proteome_UP000000807.xml.gz"
DEFAULT_TAXID=243273

from docopt import docopt
from pyproteinsExt import uniprot as pExt
from . import Unigo as createGOTreeTest
from . import uloads as createGOTreeTestFromAPI

from . import Univgo as createGOTreeTestUniverse
from . import vloads as createGOTreeTestUniverseFromAPI

from .api.store import bootstrap as goStoreStart

from .api.store.client import addTreeByTaxid as goStoreAdd
from .api.store.client import delTreeByTaxid as goStoreDel
from .api.store.client import handshake

from .stat_utils import applyOraToVector
from . import utils

from .api.pwas import listen as pwas_listen

if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)

    nDummy      = int(arguments['--size'])
    nTop        = int(arguments['--head'])
    goApiPort   = arguments['--goport'] 
    goApiHost   = arguments['--gohost'] 
    pwasApiPort = arguments['--pwasport']
    method      = arguments['--method']

    proteomeXML = arguments['--prot'] if arguments['--prot'] else DEFAULT_PROTEOME
    taxid = arguments['--taxid'] if arguments["--taxid"] else DEFAULT_TAXID
    owlFile = arguments['--onto'] if arguments['--onto'] else None

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

    if arguments['store']:
        if arguments['add'] or arguments['del']: # Client API
            handshake(goApiHost, goApiPort)
        try :
            if not arguments['dry']: # Just resume service dont alter DBB content
                    uColl = pExt.EntrySet(collectionXML=proteomeXML)
                    _ = uColl.taxids
                    if len(_) != 1:
                        raise ValueError(f"Taxids count is not equal to 1 ({len(_)}) in uniprot collection : {_}")
                    _ = _[0]

                    if arguments['del']:
                        goStoreDel([_])
                        exit(0)

                    tree_universe = createGOTreeTestUniverse( 
                                                    owlFile     = owlFile,
                                                    ns          = "biological process", 
                                                    fetchLatest = False,
                                                    uniColl     = uColl)
                except Exception as e:
                    print(f"Fatal error {e}, early Exit")
                    exit(1)

        if arguments['add']:   
            goStoreAdd(trees=[tree_universe], taxids=[_])
        elif arguments['start'] or arguments['wipe'] or arguments['dry']:
            if arguments["dry"]:
                print(f"Resuming goStore service")
                _trees  = None
                _taxids = None    
            else:
                _trees  = [tree_universe]
                _taxids = [_]
                print(f"Starting goStrore service with proteome{proteomeXML} as a universe GO tree")
            
            app = goStoreStart(trees=_trees, taxids=_taxids,\
                wipe = True if arguments['wipe'] else False,\
                cacheType='local' if not arguments['--redis'] else 'redis',\
                rp=arguments['--rp'], rh=arguments['--rh'],\
                _main_ = __name__ == '__main__' )  
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
            pwas_app = pwas_listen(goApiPort, arguments['--vectorized'])
            pwas_app.run(debug=True, port=pwasApiPort)
        else: # cli
            if not arguments['--vectorized']:# Tree ora       
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
            else: # Vector ora     
                import json
                print("Running the vectorized ora")
                resp = utils.unigo_vector_from_api(goApiPort, taxid)
                if resp.status_code != 200:
                    print(f"request returned {resp.status_code}")
                
                vectorizedProteomeTree = json.loads(resp.text)
                res = applyOraToVector(vectorizedProteomeTree, expUniprotID, deltaUniprotID, 0.5)
                print(res)
