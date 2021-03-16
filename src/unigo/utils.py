import requests

PROXIES = {
  'http': '',
  'https': '',
}

def check_proteins_subset(major_list:list[str], sublist:list[str]) -> bool:
    return set(sublist).issubset(major_list)

def unigo_tree_from_api(api_port:int, taxid:int) -> str:
    '''Interrogate GO store API and return requests response'''
    go_url = f"http://127.0.01:{api_port}/unigo/{taxid}"
    print(f"Interrogate {go_url} for go tree")
    return requests.get(go_url, proxies = PROXIES)
    
