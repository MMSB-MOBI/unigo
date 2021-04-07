import requests

HOSTNAME=None
PORT=None

def handshake(hostname, port):
    try:
        url = f"http://{hostname}:{port}/handshake"
        req = requests.get(url)
    except ConnectionError as e:
        raise ConnectionError(f"Unable to handshake at {url}\n->{e}")

    if not req.status_code == requests.codes.ok:
        raise ConnectionError(f"Error {req.status_code} while handshaking at {url}")
        if not req.text == "Hello world":
            raise ConnectionError(f"Improper handshake ({req.text}) at {url}")
    
    global HOSTNAME, PORT
    
    HOSTNAME = hostname 
    PORT     = port
    return True

def addTreeByTaxid(trees, taxids):
    print("Adding stuff")
    # Srializing it
    for tree, taxid in zip(trees, taxids):
        d = tree.serialize()
        url = f"http://{HOSTNAME}:{PORT}/add/unigo/{taxid}"
        req = requests.post(url, json=d)
        if req.status_code == requests.codes.ok:
            print(f"Successfull tree adding at {url}")
        else:
            print(f"Error {req.status_code} while inserting at {url}")

def delTreeByTaxid(taxids):
    print(f"Want to del by taxids {taxids}")
    for taxid in taxids:
        url = f"http://{HOSTNAME}:{PORT}/del/unigo/{taxid}"
        req = requests.delete(url)
        if req.status_code == requests.codes.ok:
            print(f"Successfull tree deleting at {url}")
        else:
            print(f"Error {req.status_code} while deleting at {url}")
