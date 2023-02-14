import requests
from pyproteinsext import uniprot as pExt
from . import Univgo as UniverseGO_tree
from .tree import enumNS as GoNamespaces
from uniprot_redis.store import UniprotStore
from .api.store.client import InsertionError

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

def unigo_tree_from_api(api_host:str, api_port:int, taxid:int) -> str:
	'''Interrogate GO store API and return requests response'''
	go_url = f"http://{api_host}:{api_port}/unigos/{taxid}"
	print(f"Interrogate {go_url} for go tree")
	return requests.get(go_url, proxies = PROXIES)

def unigo_vector_from_api(api_host:str, api_port:int, taxid:int) -> str:
	'''Interrogate GO store API and return requests response'''
	go_url = f"http://{api_host}:{api_port}/vectors/{taxid}"
	return requests.get(go_url, proxies = PROXIES)

def unigo_culled_from_api(api_host:str, api_port:int, taxid:int, goParameters:{}):
	'''Interrogate GO store API and return requests response'''
	go_url = f"http://{api_host}:{api_port}/vectors/{taxid}"
	print(f"Interrogate {go_url} for go culled vector {goParameters}")
	return requests.post(go_url, proxies = PROXIES, json=goParameters)

"""
	List current collection on uniprotstore

"""

def get_available_uniprot_collection(h:str,p:int):
	store = UniprotStore(host=h, port=p)

	avail_coll = store.list_collection()
	return { t[0] : len(t[1])  for t in avail_coll }
	
def sync_to_uniprot_store(host:str, port:int, owlFile:str, opt_coll_key=None):
	coll_view = get_available_uniprot_collection(host, port)
	coll_repr =  "\n".join([ f"{k}:{v} entries" for k,v in coll_view.items()] )
	user_coll_id = None
	if opt_coll_key:
		if not opt_coll_key in coll_view:
			print(coll_view)
			raise KeyError(f"{opt_coll_key} is not a collection of:\n" + coll_repr)
		user_coll_id = opt_coll_key
	else :
		print("Following uniprot collections were found")
		print(coll_repr)

		while True:
			try:
				user_coll_id = input("Please specify a collection: ")
				_ = coll_view[user_coll_id]
			except KeyError:
				print("Sorry, I didn't understand that.")
				continue
			break
	print(f"Processing Annotation trees over {user_coll_id} collection")
	#print(avail_coll)
	for ns in GoNamespaces:
		uColl = store. get_protein_collection(user_coll_id)
		print(f"\tExtracting following GO ns {ns}\n\tThis may take a while...")
		tree_universe = UniverseGO_tree( 
										owlFile     = owlFile,
										ns          = ns,#"biological process", 
										fetchLatest = False,
										uniColl     = uColl)
		yield(user_coll_id, ns, tree_universe)

# MUST REPLACE THIS W/ UNIPROT STORE
def loadUniprotContainerIter(proteomeXML, strict=True):
	"""Parse provided Uniprot XML and return following 2-uple 
		: (<first_taxid_in_collection>, <Uniprot_Collection>)
	"""
	print(f"Loading a protein collection from {proteomeXML}")

	uColl = pExt.EntrySet(collectionXML=proteomeXML)
	_ = uColl.taxids
	if len(_) != 1 and strict:
		raise ValueError(f"Taxids count is not equal to 1 ({len(_)}) in uniprot collection : {_}")
	if len(_) != 1 and not strict:
		print(f"Warning: Taxids count is not equal to 1 ({len(_)}) in uniprot collection : {_}")
	print(f"Loaded uniprot collection taxid(s){_}")

	# Patching for iterator over protein container uniprot_redis like
	def wrap():
		print("Wraped" + str(uColl))
		print(uColl.isXMLCollection)
		print(uColl.keys())
		for uniprot_id in uColl.keys():
			e = uColl.get(uniprot_id)
			d = e.toJSON()
			setattr(d, 'go', d['go'])
			print("==>" + d.go)
			yield d
	it =  wrap()
	return _[0], it#uColl

    
def generateDummySet(uColl, nTotal, deltaFrac = 0.1):
	nDummy = int(deltaFrac * nTotal)
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
	""" Yields 3 consecutive universal UniGO trees (one per GO namespace)
		foreach provided proteome XML ressource
		The iterator is a 3-uple
		(taxid:str, namespace:string, tree:Unigo.tree)
	"""
	#print(f"Loading following XML proteome(s) {proteomeXMLs}\n This may take a while...")

	uTaxids = []
	uTrees  = []
	for proteomeXML in proteomeXMLs:
		print(f"Loading following XML proteome(s) {proteomeXML}\n This may take a while...")
		uTaxid, uColl = loadUniprotContainerIter(proteomeXML)
		#uTaxids.append(uTaxid)
	
		for ns in GoNamespaces:
			print(f"\tExtracting following GO ns {ns}\n\tThis may take a while...")
			tree_universe = UniverseGO_tree( 
											owlFile     = owlFile,
											ns          = ns,#"biological process", 
											fetchLatest = False,
											uniColl     = uColl)
			yield(uTaxid, ns, tree_universe)


		#uTrees.append(tree_universe)
	
	#return zip(uTaxids, uTrees)


def parseGuessTreeIdentifiers(identifersList):
	""" Returns a list of tree key identifiers guessed from input list
		Input list elements can be actual tree keys or proteomeXML files
	"""
	guessedTreeID = []
	for ressourceMaybeXML in identifersList:
		try:
			_, uColl = loadUniprotContainerIter(ressourceMaybeXML)
			guessedTreeID += uColl.taxids
		except FileNotFoundError:
			#print(f"{ressourceMaybeXML} does not look like a file")
			guessedTreeID.append(ressourceMaybeXML)
		except IsADirectoryError:
			print(f"Warning {ressourceMaybeXML} looks like a directory")
			guessedTreeID.append(ressourceMaybeXML)
		except ValueError:
			raise IOError(f"{ressourceMaybeXML} looks like an invalid proteome XML file")
		#except Exception as e:
		#	print(type(e).__name__)
	return list(set(guessedTreeID))