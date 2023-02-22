"""Go ontology tree manipulation tool and microservice

Usage:
    unigo store cli
    unigo store server start [clear] [--rh=<redis_host> --rp=<redis_port> --gp=<store_port>]   
    unigo store client sync <owlFile> [--gp=<store_port> --gh=<store_host> --rp=<redis_port> --rh=<redis_host>] 
    unigo store client del <uniprot_coll_id> [--gp=<store_port> --gh=<store_host>]
    unigo pwas server start (vector|tree) [--pwp=<pwas_port> --gp=<store_port> --gh=<store_host>]
    unigo pwas client test <uniprot_coll_id> <n_exp> <f_delta> [--meth=<stat> --pwh=<pwas_host> --pwp=<pwas_port> --rp=<redis_port> --rh=<redis_host> --head=<n_best_pvalue>] 
    unigo pwas compute (vector|tree) <taxid> <expressed_protein_file> <delta_protein_file> [--gp=<store_port> --gh=<store_host> --method=<statMethod> --head=<n_best_pvalue>]

Options:
  -h --help     Show this screen.
  <owlFile>  ontology file location in owl format.
  --size=<number>  number of uniprotID to build the experimental protein group from the uniprot elements [default: 50]
  --delta=<number>  fraction of the experimental protein group to build the group of over/under-represented protein [default: 0.1]
  --silent  stop ORA scoring dump
  --head=<n_best_pvalue>  display n best GO pathway [default: 15]
  --gp=<store_port>  port for GO API [default: 1234]
  --gh=<host>  host name for GO API [default: localhost]
  --pwp=<pwas_port>  port for pwas API [default: 5000].
  --pwh=<pwas_host>  port for pwas API [default: localhost].
  --met=<stat>  statistical method to compute pathway p-value [default: fisher]
  <expressed_protein_file>  txt file with all proteomics experience proteins accession (one per line)
  <delta_protein_file>  txt file with all differentially expressed proteins accession (one per line)
  <uniprot_coll_id> for testing purpose, name of the proteome to use as test case generator
  <n_exp> for testing purpose, total number of observed protein
  <f_delta> for testing purpose, fraction of observed protein considered as over/under abundant
  --rp=<redis_port>  redis DB TCP port [default: 6379]
  --rh=<redis_host>  redis DB http adress [default: localhost]
"""

#unigo pwas test local  (tree|fisher|convert) [(<xmlProteomeFile> <owlFile>)] [--size=<n_exp_proteins> --delta=<n_modified_proteins> --head=<n_best_pvalue>]
    
import os

from docopt import docopt


#from . import vloads as createGOTreeUniverseFromAPI

from .api.store import bootstrap as goStoreStart

from .api.store.client.mutators import delTaxonomy as goStoreDel, addTree3NSByTaxid as goStoreAdd
from .api.store.client import handshake, get_host_param

from .utils.uniprot import sync_to_uniprot_store
from .api.pwas import listen as pwas_listen
from .repl import run as runInRepl

from .command_line import run as runSingleComputation

from .test.pwas import test_pwas_api
if __name__ == '__main__':
    arguments = docopt(__doc__)
    #print(arguments)

    goApiPort   = int(arguments['--gp'])
    goApiHost   = str(arguments['--gh'])
    redisPort   = int(arguments['--rp'])
    redisHost   = str(arguments['--rh'])
    pwasApiPort = int(arguments['--pwp'])
    method      = arguments['--method']
   
    if not arguments['server'] or ( arguments['server'] and arguments['pwas']):
        print("Configuring unigo client interface...")
        try :
            handshake(goApiHost, goApiPort)
            print(f"Successfull at {get_host_param()}")
        except Exception as e:
            print(f"Failed at {get_host_param()}")
            exit(1)

    if arguments['cli']:
       runInRepl()
    
    if arguments['store']:
        if arguments['server']: # Provider-Ressource bootstrap API
            app = goStoreStart(newElem=None,\
                clear = True if arguments['clear'] else False,\
                rp=arguments['--rp'], rh=arguments['--rh'],\
                _main_ = __name__ == '__main__' )  
            app.run(debug=False, port=goApiPort)

        elif arguments['client']: # Client-Ressource mutation API
            
            if arguments['del']:
                goStoreDel(arguments['del'])
                exit(0)
            if arguments['sync']:
                taxidTreeIter = sync_to_uniprot_store(str(redisHost), int(redisPort), str(arguments["<owlFile>"]) )
                goStoreAdd(taxidTreeIter)

    elif arguments['pwas']:
        if arguments['server']:
            pwas_app = pwas_listen(arguments['vector'])
            pwas_app.run(debug=True, port=pwasApiPort)

        elif arguments['compute']:
            print(arguments)
           
            runSingleComputation(
                arguments["<expressed_protein_file>"],\
                arguments["<delta_protein_file>"],\
                arguments["<taxid>"],\
                asVector=arguments['vector'] # Check this test
            )

        elif arguments['test']:
            print("----== Testing PWAS API ==----")
            print(f"Generating random protein sets...\n"
                  f"Reference proteome {arguments['<uniprot_coll_id>']}\n"
                  f"Protein counts [observed/modified abundance]: "
                  f"{arguments['<n_exp>']}/{arguments['<f_delta>']}"
                  )
           

            print(       arguments['--pwh'], arguments['--pwp'], \
                        arguments['<uniprot_coll_id>'],\
                        int(arguments['<n_exp>']), float(arguments['<f_delta>']),\
                        redisHost, redisPort,
                        arguments['--head'])

            test_pwas_api(arguments['--pwh'], arguments['--pwp'], \
                        arguments['<uniprot_coll_id>'],\
                        int(arguments['<n_exp>']), float(arguments['<f_delta>']),\
                        redisHost, redisPort,
                        n_top = arguments['--head'])
