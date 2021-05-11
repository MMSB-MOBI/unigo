import requests, time

HOSTNAME=None
PORT=None

class ClientError(Exception):
    def __init__(self, url):
        self.url = url

class BuildConnectionError(ClientError):
    def __init__(self, url, code):
        super().__init__(url)
        self.code = code
    def __str__(self):
        return f"Error {self.code} while monitoring build process [{self.url}]"

class InsertionError(ClientError):
    def __init__(self, url, code):
        super().__init__(url)
        self.code = code
    def __str__(self):
        return f"Insertion was denied, The trees may already exist in database [{self.url}]"

class DeletionError(ClientError):
    def __init__(self, url, code):
        super().__init__(url)
        self.code = code
    def __str__(self):
        return f"Deletion was denied, The trees may not exist in database [{self.url}]"

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

def addTree3NSByTaxid(treeTaxidIter, fromCli=False):
    
    requestedTree = {}
    for taxid, _, tree in treeTaxidIter:
        requestedTree[ f"{taxid}:{tree.ns}" ] = tree.serialize()

    url = f"http://{HOSTNAME}:{PORT}/add/taxid/{taxid}"
    req = requests.post(url, json=requestedTree)
    
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
    for taxid in taxids:
        url = f"http://{HOSTNAME}:{PORT}/del/taxid/{taxid}"
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
    url = f"http://{HOSTNAME}:{PORT}/build/vectors"
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



def unigoList(_elem="all"):
    
    d = {}
    for elem in ("trees", "vectors", "culled"):
        if _elem == "all" or elem  == _elem :
            url = f"http://{HOSTNAME}:{PORT}/list/{elem}"
            req = requests.get(url)
            d[elem] = req.json()[elem]
    return d