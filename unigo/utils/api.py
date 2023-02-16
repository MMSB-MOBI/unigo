import requests

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
	print(f"Interrogate {go_url} for go plain vector")
	return requests.get(go_url, proxies = PROXIES)

def unigo_culled_from_api(api_host:str, api_port:int, taxid:int, goParameters:{}):
	'''Interrogate GO store API and return requests response'''
	go_url = f"http://{api_host}:{api_port}/vectors/{taxid}"
	print(f"Interrogate {go_url} for go culled vector {goParameters}")
	return requests.post(go_url, proxies = PROXIES, json=goParameters)
