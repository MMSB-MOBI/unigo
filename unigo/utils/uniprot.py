from .. import Unigo
from ..tree import enumNS as GoNamespaces
from uniprot_redis.store import UniprotStore

"""
	List current collection on uniprotstore

"""

def get_view_available_uniprot_collection(h:str,p:int):
	store = UniprotStore(host=h, port=p)

	avail_coll = store.list_collection()
	d = { t[0] : len(t[1])  for t in avail_coll }
	s = "\n".join([ f"{k}:{v} entries" for k,v in d.items()] )
	return (d, s)

def sync_to_uniprot_store(host:str, port:int, owlFile:str, opt_coll_key=None):
	coll_view, coll_repr = get_view_available_uniprot_collection(host, port)
	
	user_coll_id = None
	if opt_coll_key:
		if not opt_coll_key in coll_view:
			print(coll_view)
			raise KeyError(f"{opt_coll_key} is not a collection of:\n" + coll_repr)
		user_coll_id = opt_coll_key
	else :
		print("Please choose one of following uniprot collections to build GO terms tree")
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
	store = UniprotStore(host=host, port=port)
	for ns in GoNamespaces:
		uColl = store.get_protein_collection(user_coll_id)
		print(f"\tExtracting following GO ns {ns}\n\tThis may take a while...")
		unigo_blueprint = Unigo( 
										owlFile     = owlFile,
										ns          = ns,#"biological process", 
										fetchLatest = False,
										uniColl     = uColl)
		yield(user_coll_id, ns, unigo_blueprint)

"""
Generate two list of uniprot identifiers __the observed and the delta abundant proteins__ from
a uniprot collection __the proteome__ specified by a collection identifier
"""
def generate_dummy_sets(coll_id, n_total, n_delta_frac, h:str,p:int, seed:int):
	coll_view, coll_repr = get_view_available_uniprot_collection(h, p)
	if not coll_id in coll_view:
		raise KeyError(f"The uniprot collection identifer you provided {coll_id} is not valid:\n", coll_repr)
	print(coll_view)
	if n_total > coll_view[coll_id]:
		raise ValueError("Specified observed protein count ({n_total}) exceeds total collection size {coll_view[coll_id]}")
	
	obs_uniprot_ids =[]
	store = UniprotStore(host=h, port=p)
	uColl = store.get_protein_collection(coll_id)
	for uObj in uColl:
		if len(obs_uniprot_ids) == n_total:
			break
		if uObj.go:
			obs_uniprot_ids.append(uObj.id)
	print(f"Managed to pull a obs total count of {len(obs_uniprot_ids)} uniprot identifiers w/ GO annotation\n")
	
	n_total = len(obs_uniprot_ids)
	n_dummy_delta = int(float(n_delta_frac) * n_total)
	print(f"{coll_id} proteome generated dummy observed/delta abundant data sets of sizes {n_total}/{n_dummy_delta}")

	return obs_uniprot_ids, obs_uniprot_ids[:n_dummy_delta]