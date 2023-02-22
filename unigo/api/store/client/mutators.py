import requests, time
from . import get_host_param
from . import InsertionError, DeletionError, BuildConnectionError

def addTree3NSByTaxid(treeTaxidIter, fromCli=False):
    hostname, port = get_host_param()
    requestedTree = {}
    for taxid, _, tree in treeTaxidIter:     
        if taxid not in requestedTree : requestedTree[taxid] = {}  
        requestedTree[taxid][ f"{taxid}:{tree.ns}" ] = tree.serialize()
    
    for taxid in requestedTree:
        url = f"http://{hostname}:{port}/add/taxid/{taxid}"
        req = requests.post(url, json=requestedTree[taxid])
        if req.status_code == requests.codes.ok:
            msg = f"Successfull tree adding at {url}"
            if fromCli:
                return msg
            print(msg)
            
        else:
            if fromCli:
                raise InsertionError(url, req.status_code)
            print(f"Error {req.status_code} while inserting at {url}") 

def delTaxonomy(taxids, fromCli=False):
    #print(f"Want to del by taxids {taxids}")
    hostname, port = get_host_param()
    for taxid in taxids:
        url = f"http://{hostname}:{port}/del/taxid/{taxid}"
        req = requests.delete(url)

        if req.status_code == requests.codes.ok:
            msg = f"Successfully deleted data under taxonomy {taxid} [{url}]"
            if fromCli:
                return msg
            print(msg)
                
        else:
            if fromCli:
                raise DeletionError(url, req.status_code)            
            print(f"Error {req.status_code} while deleting at {url}")


def buildVectors(fromCli=False):
    """ Trigger Vector building
        * List total number of vector to build
        * regularly asks for status
    """
    (status, size) = _pingAndUnwrapBuildReq(fromCli)
    yield (status, size)
    while status in ["starting", "running"]:
        time.sleep(1)
        (status, _size) = _pingAndUnwrapBuildReq(fromCli)
        if _size < size:
            for iSize in range(size, _size, -1):
                time.sleep(0.5)
                if iSize - 1 > 0:
                    yield (status, iSize - 1)
            size = _size
    yield ("completed", 0)


def _pingAndUnwrapBuildReq(fromCli):
    hostname, port = get_host_param()
    url = f"http://{hostname}:{port}/build/vectors"
    req = requests.get(url)
    if not req.status_code in ["200", "202"]:
        data = req.json()
        status = data["status"]
        n = 0 if status == "nothing to build" else len(data["targets"])
        return (status, n)
    else:
        if fromCli:
            raise BuildConnectionError(url, req.status_code)            
        print(f"Error {req.status_code} while buidling at {url}")
