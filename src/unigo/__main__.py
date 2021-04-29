"""Go ontology tree manipulation tool and microservice

Usage:
    unigo store server redis start [clear] [<owlFile> <xmlProteomeFile>... ] [--rh=<redis_host> --rp=<redis_port> --go=<store_port>]
    unigo store server local start <xmlProteomeFile>... [--gp=<store_port>]
    unigo store client add  <owlFile> <xmlProteomeFile>... [--gp=<store_port> --gh=<host>]
    unigo store client del <xmlProteomeFile>... [--gp=<portNumber> --gh=<hostname>]  
    unigo pwas server (vector|tree) [--pwp=<pwas_port> --gp=<store_port> --gh=<store_host>]
    unigo pwas cli (vector|tree) <taxid> <expressed_protein_file> <delta_protein_file> [--gp=<store_port> --gh=<store_host> --method=<statMethod>]
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
  --pwp=<pwas_port>  port for pwas API [default: 5000].
  --method=<statMethod>  statistical method to compute pathway p-value [default: fisher]
  --exp-prot=<list.txt>  txt file with all proteomics experience proteins accession (one per line)
  --de-prot=<list.txt>  txt file with all differentially expressed proteins accession (one per line)
  --vectorized  use the goStore vectorize protocol
  --redis  use redis backend as storage engine
  --rp=<redis_port>  redis DB TCP port [default: 6379]
  --rh=<redis_host>  redis DB http adress [default: localhost]
"""
import os

from docopt import docopt


#from . import vloads as createGOTreeUniverseFromAPI

from .api.store import bootstrap as goStoreStart

from .api.store.client import addTreeByTaxid as goStoreAdd
from .api.store.client import delTreeByTaxid as goStoreDel
from .api.store.client import handshake


from .utils import loadUniversalTreesFromXML

from .api.pwas import listen as pwas_listen

if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)

    nDummy      = int(arguments['--size'])
    nTop        = int(arguments['--head'])
    goApiPort   = arguments['--gp'] 
    goApiHost   = arguments['--gh'] 
    pwasApiPort = arguments['--pwp']
    method      = arguments['--method']
       
    #taxid = arguments['--taxid'] if arguments["--taxid"] else DEFAULT_TAXID
    #owlFile = arguments['--onto'] if arguments['--onto'] else None

 #   if arguments['fisher'] or arguments['convert']:
 #       arguments['tree'] = True



    #Load or create proteins sets
    # With default inputs for proteome and annotations
    #              settings 
    # TO RUN & VALIDATE
    if arguments['test']:
        target = list(set(['fisher', 'convert', 'tree']) & set(arguments.keys()))[0]
        test.run(nDummy, nTop, target=target)
        exit(0)
   
    
        
    

    if arguments['server']: # Provider-Ressource bootstrap API
        if not arguments["<xmlProteomeFile>"]:
            print(f"Resuming goStore service")
            taxidTreeIter = None 
        else:
            print(f"Starting goStore service with proteome{arguments['<xmlProteomeFile>']} as a universe GO tree")
            taxidTreeIter = loadUniversalTreesFromXML(\
                arguments["<xmlProteomeFile>"],\
                arguments["<owlFile>"])
        
        app = goStoreStart(newElem=taxidTreeIter,\
            clear = True if arguments['clear'] else False,\
            cacheType='local' if not arguments['redis'] else 'redis',\
            rp=arguments['--rp'], rh=arguments['--rh'],\
            _main_ = __name__ == '__main__' )  
        app.run(debug=False, port=goApiPort)


    if arguments['client']: # Client-Ressource mutation API
        handshake(goApiHost, goApiPort)
        if arguments['del']:
            goStoreDel([_])
            exit(0)
        if arguments['add']:
            goStoreAdd(trees=[tree_universe], taxids=[uTaxid])
    

                    
        #except Exception as e:
        #    print(f"Fatal error {e}, early Exit")
        #    exit(1)

    if arguments['pwas']:
        if arguments['server']:
            pwas_app = pwas_listen(goApiPort, arguments['vector'])
            pwas_app.run(debug=True, port=pwasApiPort)

        elif arguments['cli']:
            cli.run(guments["expressed_protein_file"],\
                    arguments["delta_protein_file"],\
                    goApiPort, taxid,
                    asVector=arguments['vector'] # Check this test
            )