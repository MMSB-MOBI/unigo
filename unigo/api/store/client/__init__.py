import os
import requests

HOSTNAME = None
PORT     = None

def get_host_param():
    if not HOSTNAME or not PORT:
        raise ConnectionError(f"You must handshake unigo server")
    return (HOSTNAME, PORT)

def handshake(hostname, port):
    global HOSTNAME, PORT

    try:
        url = f"http://{hostname}:{port}/ping"
        req = requests.get(url)
    except ConnectionError as e:
        raise ConnectionError(f"Unable to handshake at {url}\n->{e}")

    if not req.status_code == requests.codes.ok:
        raise ConnectionError(f"Error {req.status_code} while handshaking at {url}")
    if not req.text == "pong":
        raise ConnectionError(f"Improper handshake ({req.text}) at {url}")
    
    HOSTNAME = hostname 
    PORT     = port
    return True


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




