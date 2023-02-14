"""Go ontology tree manipulation tool and microservice

Usage:
    unigo store cli
    unigo store server start [clear] [--rh=<redis_host> --rp=<redis_port> --go=<store_port>]   
    unigo store client sync <owlFile> [--coll=<protein_collections> --gp=<store_port> --gh=<store_host> --up=<uniprot_redis_port> --uh=<uniprot_redis_host>] 
    unigo store client del <coll_id> [--gp=<store_port> --gh=<store_host>]
    unigo pwas server (vector|tree) [--pwp=<pwas_port> --gp=<store_port> --gh=<store_host>]
    unigo pwas compute (vector|tree) <taxid> <expressed_protein_file> <delta_protein_file> [--gp=<store_port> --gh=<store_host> --method=<statMethod>]
    unigo pwas test local  (tree|fisher|convert) [(<xmlProteomeFile> <owlFile>)] [--size=<n_exp_proteins> --delta=<n_modified_proteins> --head=<n_best_pvalue>]
    unigo pwas test client [--port=<portNumber>] [--method=<statMethod>] [--prot=<xmlFile>] 

Options:
  -h --help     Show this screen.
  <owlFile>  ontology file location in owl format.
  <xmlProteomeFile>  uniprot file location in xml format.
  --size=<number>  number of uniprotID to build the experimental protein group from the uniprot elements [default: 50]
  --delta=<number>  fraction of the experimental protein group to build the group of over/undee-represented protein [default: 0.1]
  --silent  stop ORA scoring dump
  --head=<n_best_pvalue>  display n best GO pathway [default: 5]
  --gp=<store_port>  port for GO API [default: 1234]
  --gh=<host>  host name for GO API [default: localhost]
  --up=<store_port>  port for Uniprot API [default: 6379]
  --uh=<host>  host name for Uniprot API [default: localhost]
  --coll=<protein_collections> protein collection identifiers separated by comma
  --pwp=<pwas_port>  port for pwas API [default: 5000].
  --method=<statMethod>  statistical method to compute pathway p-value [default: fisher]
  <expressed_protein_file>  txt file with all proteomics experience proteins accession (one per line)
  <delta_protein_file>  txt file with all differentially expressed proteins accession (one per line)
  --rp=<redis_port>  redis DB TCP port [default: 6379]
  --rh=<redis_host>  redis DB http adress [default: localhost]
"""
import os

from docopt import docopt


#from . import vloads as createGOTreeUniverseFromAPI

from .api.store import bootstrap as goStoreStart

from .api.store.client import addTree3NSByTaxid as goStoreAdd
from .api.store.client import delTaxonomy as goStoreDel
from .api.store.client import handshake

from .utils import sync_to_uniprot_store
from .api.pwas import listen as pwas_listen
from .repl import run as runInRepl

from .command_line import run as runSingleComputation

if __name__ == '__main__':
    arguments = docopt(__doc__)
    #print(arguments)

    nDummy      = int(arguments['--size'])
    nTop        = int(arguments['--head'])
    goApiPort   = arguments['--gp'] 
    goApiHost   = arguments['--gh'] 
    uniprotApiPort   = arguments['--up'] 
    uniprotApiHost   = arguments['--uh'] 
    pwasApiPort = arguments['--pwp']
    method      = arguments['--method']

    if arguments['cli']:
       runInRepl()
    
    if arguments['store']:
        #Load or create proteins sets
        # With default inputs for proteome and annotations
        #              settings 
        # TO RUN & VALIDATE
        if arguments['test']:
            target = list(set(['fisher', 'convert', 'tree']) & set(arguments.keys()))[0]
            test.run(nDummy, nTop, target=target)
            exit(0)
    
        if arguments['server']: # Provider-Ressource bootstrap API
            app = goStoreStart(newElem=None,\
                clear = True if arguments['clear'] else False,\
                cacheType='local' if not arguments['redis'] else 'redis',\
                rp=arguments['--rp'], rh=arguments['--rh'],\
                _main_ = __name__ == '__main__' )  
            app.run(debug=False, port=goApiPort)

        elif arguments['client']: # Client-Ressource mutation API
            handshake(goApiHost, goApiPort)
            if arguments['del']:
                goStoreDel(arguments['del'])
                exit(0)
            if arguments['sync']:
                taxidTreeIter = sync_to_uniprot_store(str(uniprotApiHost), int(uniprotApiPort), str(arguments["<owlFile>"]) )
                goStoreAdd(taxidTreeIter)

    elif arguments['pwas']:
        handshake(goApiHost, goApiPort)
        if arguments['server']:
            pwas_app = pwas_listen(goApiHost, goApiPort, arguments['vector'])
            pwas_app.run(debug=True, port=pwasApiPort)

        elif arguments['compute']:
            print(arguments)
           
            runSingleComputation(
                arguments["<expressed_protein_file>"],\
                arguments["<delta_protein_file>"],\
                goApiHost,\
                goApiPort,\
                arguments["<taxid>"],\
                asVector=arguments['vector'] # Check this test
            )