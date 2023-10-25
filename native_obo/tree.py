import networkx as nx
from typing import Iterator, Union, NewType, Literal, get_args
from networkx.classes.digraph import DiGraph
from .io_obo import obo_node_buffer_iter
from uniprot_redis.store.schemas import UniprotDatum, UniprotAC
from pyproteinsext.uniprot import Entry as Uniprot

from decorator import decorator

NodeID = NewType("NodeID", str)
ProteinsType = Literal["background", "target"]

# Move this to decorator
def literal_assert(v:str, lit:Literal)->None:
    assert v in get_args(lit), (f"Value {v} is not part of literal {get_args(lit)}")

# Somehow painfulls
#@decorator
#def literal_arg_checker(f, *args, **kwargs):
#    for 

class GO_tree(nx.DiGraph):
    def __init__(self):
        super().__init__()
        self.protein_load_status = (False, False)

    def add_node(self, id, **params):
        """
        nx add_node wrapper to account for node aliases
        """
        super(GO_tree, self).add_node(id, **params)
        if 'alt_id' in params:
            alt_ids = params['alt_id'] if type(params['alt_id']) is list else [ params['alt_id'] ]        
            real_n = self.nodes[id]
            for alt_id in alt_ids:
                super(GO_tree, self).add_node(alt_id, alias_to = real_n, _id = alt_id)
    
    def get_go_node(self, node_id:NodeID):
        """
        nx G.node wrapper for the automatic forward of alias go_id to concrete node 
        """
        _ = self.nodes[node_id]
        if "alias_to" in _:
            return _["alias_to"]
        return _

    def concrete_nodes(self)->Iterator[dict]:
        for k in self.nodes:
            d = self.maybe_concrete(k)
            if d is None:
                continue
            yield d

    def maybe_concrete(self, node_id:NodeID)->Union[dict, None]:
        """ Boolean evaluation of queried node as a node that is not obsolete nor and alias
        """
        _ = self.nodes[node_id]
        if not 'is_obsolete' in _ and not 'alias_to' in _:
            return _
        return None
        
    def successors(self, node_id:NodeID):
        """ Wrapping for alias forwarding 
            In case we wanna ascend from a deprecated GO term
        """
        _ = self.nodes[node_id]
        if "alias_to" in _:
            return super(GO_tree, self).successors(_["alias_to"]["_id"])
        return super(GO_tree, self).successors(node_id)

    def load_proteins(self, k:ProteinsType, protein_coll:Iterator[ Union[UniprotDatum, Uniprot] ]) -> None:
        """
        Iterate through a uniprot datum collection and attach UniprotDatum to 
        corresponding concerte node to background or target attribute
        """
        literal_assert(k, ProteinsType)
        self.clear_proteins(k)
        for unip_datum in protein_coll:
            for go_datum in unip_datum.go:
                try:
                    go_node = self.get_go_node(go_datum.id)
                    if not k in go_node:
                        go_node[k] = set()
                    go_node[k].add(unip_datum)
                except KeyError:
                    print (f"{go_datum.id} not found")
        self.protein_load_status = (True, self.protein_load_status[1]) if k == "background" \
        else (self.protein_load_status[0], True)
 
    def clear_proteins(self, k: ProteinsType):
        literal_assert(k, ProteinsType)
        for node_dic in self.nodes.values():
            _ = node_dic.pop(k, None)
        self.protein_load_status = (False, self.protein_load_status[1]) if k == "background" \
        else (self.protein_load_status[0], False)

    @property
    def root_ids(self)->set[NodeID]:
        r = set()
        for k in self.nodes:
            d = self.get_go_node(k)
            if not '_is_a' in d and not 'is_obsolete' in d:
                r.add(d['_id'])
        return r

    def get_proteins(self, node_id:NodeID, k:ProteinsType="background", deep=True)->set[UniprotDatum]:
        """ Get all UniprotDatum attached to subtree rooted at provided node id """
        literal_assert(k, ProteinsType)
        def _get_proteins(node_id, deep)->set[UniprotDatum]:
            curr_node = self.get_go_node(node_id)
            results   = curr_node[k] if k in curr_node else set()
            if not deep:
                return results
            
            for next_node_id in self.successors(node_id):
                results |= _get_proteins(next_node_id, deep)
            return results
            
        return _get_proteins(node_id, deep)
    
    @property
    def leave_ids(self)->Iterator[NodeID]:
        for d in self.concrete_nodes():
            bot = True
            for _ in self.successors(d['_id']):
                bot = False
                break
            if bot:
                yield d['_id']
    
    @property
    def ora_rdy(self):
        return self.protein_load_status[0] and self.protein_load_status[1] 

    #def cut(self: node:NodeID)->Self:
    # mmmh

# MAybe move to io
def reader(obo_file_path:str, keep_obsolete=True, ns=None)-> DiGraph :
    G = GO_tree() #nx.DiGraph()

    for node_buffer in obo_node_buffer_iter(open(obo_file_path, 'r')):
        #print(node_buffer)
        #print(node_buffer['id'])
        #print(node_buffer.nx_node_param)
        G.add_node(node_buffer['id'], **node_buffer.nx_node_param)
        for node_parent_id in node_buffer.is_a_iter():
            G.add_edge(node_parent_id, node_buffer['id'], type='is_a')

        if not keep_obsolete or not  node_buffer.is_obsolete:
            continue
        for ref_node_id in node_buffer.consider_iter():
            G.add_edge(node_buffer['id'], ref_node_id, type='refer_to')

        #print(G.nodes[ node_buffer['id'] ])
        #print(f"==>{G[ node_parent_id ][]}")
        
        
    return G

