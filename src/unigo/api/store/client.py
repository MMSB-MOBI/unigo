import requests

HOSTNAME=None
PORT=None

def handshake(hostname, port):
    try:
        url = f"http://{hostname}:{port}/ping"
        req = requests.get(url)
    except ConnectionError as e:
        raise ConnectionError(f"Unable to handshake at {url}\n->{e}")

    if not req.status_code == requests.codes.ok:
        raise ConnectionError(f"Error {req.status_code} while handshaking at {url}")
        if not req.text == "pong":
            raise ConnectionError(f"Improper handshake ({req.text}) at {url}")
    
    global HOSTNAME, PORT
    
    HOSTNAME = hostname 
    PORT     = port
    return True

def addTree3NSByTaxid(treeTaxidIter):
    
    requestedTree = {}
    for taxid, _, tree in treeTaxidIter:
        requestedTree[ f"{taxid}:{tree.ns}" ] = tree.serialize()

    url = f"http://{HOSTNAME}:{PORT}/add/taxid/{taxid}"
    req = requests.post(url, json=requestedTree)
    
    if req.status_code == requests.codes.ok:
        print(f"Successfull tree adding at {url}")
    else:
        print(f"Error {req.status_code} while inserting at {url}")

def delTaxonomy(taxids):
    print(f"Want to del by taxids {taxids}")
    for taxid in taxids:
        url = f"http://{HOSTNAME}:{PORT}/del/taxid/{taxid}"
        req = requests.delete(url)
        if req.status_code == requests.codes.ok:
            print(f"Successfully deleted data under taxonomy {taxid}[{url}]")
        else:
            print(f"Error {req.status_code} while deleting at {url}")


def buildVectors():
    """ Trigger Vector building
        * List total number of vector to build
        * regularly asks for status
    """
    url = f"http://{HOSTNAME}:{PORT}/list/trees"
    req = requests.delete(url)

def unigoList(_elem="all"):
    
    d = {}
    for elem in ("trees", "vectors", "culled"):
        if _elem == "all" or elem  == _elem :
            url = f"http://{HOSTNAME}:{PORT}/list/{elem}"
            req = requests.get(url)
            d[elem] = req.json()[elem]
    return d