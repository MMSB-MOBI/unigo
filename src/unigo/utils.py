import requests
from pyproteinsExt import uniprot as pExt
from . import Univgo as createGOTreeUniverse

PROXIES = {
	'http': '',
	'https': '',
}

def setProxy(http, https):
	global PROXIES
	PROXIES["http"] = http
	PROXIES["https"] = https
  
def check_proteins_subset(major_list:list[str], sublist:list[str]) -> bool:
	return set(sublist).issubset(major_list)

def unigo_tree_from_api(api_port:int, taxid:int) -> str:
	'''Interrogate GO store API and return requests response'''
	go_url = f"http://127.0.01:{api_port}/unigo/{taxid}"
	print(f"Interrogate {go_url} for go tree")
	return requests.get(go_url, proxies = PROXIES)
    

def unigo_vector_from_api(api_port:int, taxid:int) -> str:
	'''Interrogate GO store API and return requests response'''
	go_url = f"http://127.0.01:{api_port}/vector/{taxid}"
	print(f"Interrogate {go_url} for go vector")
	return requests.get(go_url, proxies = PROXIES)


def loadUniprotCollection(proteomeXML, strict=True):
	print(f"Loading a protein collection from {proteomeXML}")

	uColl = pExt.EntrySet(collectionXML=proteomeXML)
	_ = uColl.taxids
	if len(_) != 1 and strict:
		raise ValueError(f"Taxids count is not equal to 1 ({len(_)}) in uniprot collection : {_}")
	if len(_) != 1 and not strict:
		print(f"Warning: Taxids count is not equal to 1 ({len(_)}) in uniprot collection : {_}")
	
	return _[0], uColl

    
def generateDummySet(uColl, nTotal, deltaFrac = 0.1):
	print(f"Setting up a dummy experimental collection of {nDummy} elements")
	expUniprotID =[]
	for uObj in uColl:
		if len(expUniprotID) == nTotal:
			break
		if uObj.isGOannot and uObj.taxid == uTaxid:
			expUniprotID.append(uObj.id)
	
	nDelta = int(nDummy*deltaFrac)
	print(f"Considering {nDelta} proteins among {nTotal} experimental as of significantly modified quantities")
	return expUniprotID[:nDelta]

def loadUniprotIDsFromCliFiles(expressedProtIDFile, deltaProtIDFile):
	expUniprotID   = []
	deltaUniprotID = []
	with open(expressedProtIDFile) as f:
		for l in f:
			expUniprotID.append(l.rstrip())
	with open(deltaProtIDFile) as f:
		for l in f:
			deltaUniprotID.append(l.rstrip())

	if not check_proteins_subset(expUniprotID, deltaUniprotID):
		raise Exception("Differentially expressed proteins are not completely included in total proteins")

	return expUniprotID, deltaUniprotID


def loadUniversalTreesFromXML(proteomeXMLs, owlFile):
	print(f"Loading following XML proteome(s) {proteomeXMLs}\n This may take a while...")

	uTaxids = []
	uTrees  = []
	for proteomeXML in proteomeXMLs:
		uTaxid, uColl = loadUniprotCollection(proteomeXML)
		uTaxids.append(uTaxid)
		tree_universe = createGOTreeUniverse( 
											owlFile     = owlFile,
											ns          = "biological process", 
											fetchLatest = False,
											uniColl     = uColl)
		uTrees.append(tree_universe)
	
	return zip(uTaxids, uTrees)