from .node import BaseNode
from .heap import CoreHeap as kNodes
from ..api.data_objects import CulledGoParametersSchema as goParameterValidator
import json
import sys
sys.setrecursionlimit(3000)

enrichment_ref_pop = ['observed', 'full_proteome']
#MAX_DEPTH=9999
CURR_SHORT_PATH=9999 
class OraNode(BaseNode):
    def __init__(self, ID, name, uniprot_ids):
        super().__init__(ID)
        self.name = name
        self.is_valid = False
        self.uniprot_ids = uniprot_ids
        #self.visited     = 0 
   
    def __str__(self):
        return  f"NODE: {self.name} {self.ID}\n" + \
                f">is_a     : {[ _.name +'|'+_.ID for _ in self.is_a ]}\n" + \
                f">children : {[ _.name +'|'+_.ID for _ in self.children ]}\n"

    """    
    def bfs(self, w_node, depth):
        global CURR_SHORT_PATH
        depth += 1
        if depth >= CURR_SHORT_PATH:
            return None
        self.depth = depth if depth < self.depth else self.depth # update distance from start
        if self == w_node:
            CURR_SHORT_PATH = self.depth
            return self.depth
        
        for n in self.children:
            _ = n.bfs(w_node, depth)
            if not _ is None:
                self.depth = _
    """
    #def lineage(self,)
    # roll up a node lineage, validating nodes along the way eventually
    # droping consecutive non-valid leaves

    def rollup(self, v_nodes, main_heap, GO_constraints:goParameterValidator, ref_pop, is_leaf_streak:bool):
       
        # the current node has other child than the one calling/rolling up
        # it can't be removed
        # but was already visited, so no need to propagate upstairs
        if self in v_nodes and self.children:
            return
        
        if self.name == 'ora_root':
            return
        
        # The current node was not previously assesed
        if not self in v_nodes:
            self.is_valid = self.validate(GO_constraints, ref_pop)
        
        v_nodes.add(self)
        
        # Let's remove a dangling and invalid node 
        if not self.children and not self.is_valid:
            #print(f"\nDeleting terminal a non valid ora node")
            #print(self)
           
            main_heap.remove(self)
            #print("<Rolling up>")
            # unregister current node from its parents 
            for p_node in self.is_a:
                p_node.children = list(set(p_node.children) - set([self])) # Could be faster if postproces by taking union any_node.children | tree.node_set
                 
            # Keep on rolling up
                p_node.rollup(v_nodes, main_heap, GO_constraints, ref_pop, is_leaf_streak)

    # Apply following constraints 
    #{'minCount': 0, 'maxFreq': 0.05, 'maxCount': 50}
    # couple of values could be considered for count and bk freq and pvalue
    def validate(self, go_param:goParameterValidator, ref_pop:str):
        #print(f"Validating {self.name}")

        if self.features['fisher'] is None:            
            return False

        count_exp = len(self.uniprot_ids)
        count_sa =  len(self.features['fisher']['xp_hits'])

        bkfq_deep = self.features['fisher']['pathway_freq_deep']
        bkfq_obs  = self.features['fisher']['pathway_freq_obs']

     

        p_deep    = self.features['fisher']['pvalue_deep']
        p_annot   = self.features['fisher']['pvalue_annot']
        pvalue    = p_deep if ref_pop == 'full_proteome' else p_annot

        if count_sa < go_param['minCount']:
            return False
        if count_sa > go_param['maxCount']:
            return False
        if  bkfq_deep > go_param['maxFreq']:
            return False
        if pvalue > go_param['pvalue']:
            return False
        
        return True
    
    # serializing pvalue_deep
    def dictify(self): # we try with no heap --> redundant nodes maybe 
        d = { 
            'ID' : self.ID,
            'name' : self.name,
            'is_valid' : self.is_valid,
            'children' : []
        }
        if self.is_valid:
            d['pvalue'] = self.features['fisher']['pvalue_deep']
            d['n_xp'] = len(self.features['fisher']['xp_obs'])
            d['n_sa'] = len(self.features['fisher']['xp_hits'])
        

        for c in self.children:
            d['children'].append( c.dictify() )
        #d['nb_children'] = len(d['children'])
        return d
    
    # update is_a or children attribute
    def update(self, type, value):
        assert type in ['is_a', 'children']
        if isinstance(value, list):
            value = set(value)
        if isinstance(value, OraNode):
            value = set([value])
        if not isinstance(value, set):
            raise ValueError(f"Cant update is_a/children attribute w/ this => {value}")
        
        if type == 'is_a':
            self.is_a = list( set(self.is_a) | value  )
            return self.is_a
        
        self.children = list( set(self.children) | value  )
        return self.children
    # Not sure how to deal with second+ encounter of a node
    def collapse(self, v_node, main_heap):
        if self in v_node:
            return self
        v_node.add(self)
        #print(f"Testing {self.name} {self.is_valid} {len(self.children)}")
        if not self.children:
            return self
        if len(self.children) == 1 and not self.is_valid:
            main_heap.remove(self)
            return self.children[0].collapse(v_node, main_heap)
        self.children = [ n.collapse(v_node, main_heap) for n in self.children ]
        return self
    
 

class EmptyOraDescent(Exception):
    pass

class EmptyOraTree(Exception):
    pass

class OraTreeAccessError(Exception):
    pass

class OraTreeWireError(Exception):
    pass

class OraTree():
    def __init__(self):
        self.root = OraNode('0000', 'ora_root', None)
        self.node_set = kNodes()    
        self.leaf_set = set()
        self.depth_marked = False

    def mark_depth(self):
        #print("OraTree:mak_depth:: Indexing depth marks...")
        for node in self.node_set:
            node.features['depth'] = False
        self.root.children[0].mark_depth(0)
        self.depth_marked = True
        #print("OraTree:mak_depth:: Done")

    def bfs(self, ID_or_node):
        
        if not self.depth_marked:
            self.mark_depth()

        if isinstance(ID_or_node, str): 
            ID_or_node = self.get_node(ID_or_node)
        elif not isinstance(ID_or_node, OraNode): 
            raise TypeError(f"Ora_tree.bfs:: provided element is not a node or a string : \"{type(ID_or_node)}\"")
        #print(f"Ora_tree.bfs:: {ID_or_node.name} bfs line::\n")
        lin = [ ID_or_node ]
        while lin[-1].features['depth'] != 1:
            _ = sorted(lin[-1].is_a, key=lambda node:node.features['depth'])
            #print(f"{lin[-1].name}, must choose in:: {[ (n.name, n.features['depth']) for n in _]}")
            lin.append(_[0])
        
        return lin
    
    def get_node(self, ID)->OraNode:
        for n in self.node_set:
            if n.ID == ID:
                return n
        raise OraTreeAccessError(f"No such key {ID}")
    
    def _collapse(self)->int:
        if not self.root.children:
            raise EmptyOraTree()
           
        n = len(self.node_set)
        visited_nodes = set()
        ns_like_node = self.root.children[0]
        ns_like_node.children = [ n.collapse(visited_nodes, self.node_set) for n in ns_like_node.children ]
        return n - len(self.node_set)
    
    # roll up all leaves lineage, validating nodes along the way eventually
    # droping consecutive/terminal non-valid leaves
    # When all leaves are ro
    def prune(self, GO_constraints:goParameterValidator, ref_pop, collapse=True):
        global enrichment_ref_pop
        assert(ref_pop in enrichment_ref_pop)

        if not self.root.children:
            raise EmptyOraTree()

        if not self.leaf_set:
            print("ora_tree:prune :: not leaf set found building it...")
            for ora_node in self.root.walk():
                if not ora_node.children:
                    self.leaf_set.add(ora_node)

        print(f"ora_tree:prune :: {self.dimensions} (nodes, leaves) based on reference pop \"{ref_pop}\"")

        ndel_by_rollup = len(self.node_set)
        visited_nodes = set()
        for lf_node in self.leaf_set:
            #print(f"Ora_tree::prune: Nibbling a new leaf :: {lf_node.name}")
            lf_node.rollup(visited_nodes, self.node_set, GO_constraints, ref_pop, True)        
        ndel_by_rollup -= len(self.node_set)
        print(f"Ora_tree:prune: {ndel_by_rollup} removed by rolling up")
        print("Nb visited / Nb left::", len(visited_nodes), len(self.node_set))
        print("Nb left but never visited::", len(self.node_set.asSet - visited_nodes))
        #if collapse:
        #   ndel_by_collapse = self._collapse()
        #   print(f"Ora_tree:prune: {ndel_by_collapse} removed by collapsing")
        
        self.leaf_set = self.leaf_set & self.node_set.asSet
        ndel_by_collapse = self.build_valid_minimal()
        return (ndel_by_rollup, ndel_by_collapse)
    
    # All valid nodes will be represented in the minimal tree
    # Their path to the root is always the shortest according 
    def build_valid_minimal(self):
        print(f"OraTree::build_valid_minimal input dimensions {self.dimensions}") 
        # mark the tree
        self.mark_depth()
        
        for n in self.node_set:
            n.features['_memo']= {
                '_children' : [],
                '_is_a'     : [],
                'picked'    : False
            }
        # Build leaves to root minimal tree    
        print(f"Building minimal tree for leaf_set {len(self.leaf_set)} elements")
        for leaf in self.leaf_set:
            
            leaf_bfs = self.bfs(leaf)

            assert(not leaf_bfs[0].features['_memo']['picked'])
            #print("##[1]", leaf_bfs[1].name)
            #print("child::", str([ _.name for _ in leaf_bfs[1].children ]))
            #print(leaf_bfs[1].features)
            #print("##[0]", leaf_bfs[0].name)
            #print("is_a::", str([ _.name for _ in leaf_bfs[0].is_a ]))
            #print(leaf_bfs[0].features)
            
            leaf_bfs[0].features['_memo']['picked'] = True
            leaf_bfs[1].features['_memo']['picked'] = True
            
            leaf_bfs[0].features['_memo']['_is_a'] = [leaf_bfs[1]]
            leaf_bfs[1].features['_memo']['_children'].append(leaf_bfs[0])

            for i, n in enumerate(leaf_bfs[1:]):
                memo_child = leaf_bfs[ i-1 ].features['_memo']
                memo_parent  = leaf_bfs[ i ].features['_memo']
                assert(len(memo_child['_is_a']) == 0)
                memo_parent['_children'].append(leaf_bfs[i-1])
                memo_child['_is_a'] = [ leaf_bfs[i] ] #.append(leaf_bfs[i])
                if memo_parent['picked']:
                    return
                memo_parent['picked'] = True

        # Add to this minimal tree, the minimal tree of validated intermediary nodes not already added by leaves to root minimal tree
        node_to_add = set([ n for n in self.node_set if n.children and n.is_valid and not n.features['_memo']['picked'] ]) 
        print(f"Augmenting the minimal tree for reminder nodes {len(self.leaf_set)} elements")
        
        for rem_node in node_to_add:
            # As a reminder node, it may be consider as a leave in the new graph
            # consecutive reminder node may be son of th
            rem_bfs = self.bfs(rem_node)

            assert(not rem_bfs[0].features['_memo']['picked'])
            
            rem_bfs[0].features['_memo']['picked'] = True
            rem_bfs[1].features['_memo']['picked'] = True
            
            rem_bfs[0].features['_memo']['_is_a'] = [rem_bfs[1]]
            rem_bfs[1].features['_memo']['_children'].append(rem_bfs[0])

            for i, n in enumerate(rem_bfs[1:]):
                assert(len(rem_bfs[ i-1 ]['_memo']['_is_a']) == 0)
                rem_bfs[ i ]['_memo']['_children'].append(rem_bfs[i-1])
                rem_bfs[ i-1 ]['_memo']['_is_a'] = [ rem_bfs[ i ] ]#.append(leaf_bfs[i])
                if rem_bfs[ i ]['_memo']['picked']:
                    return
                rem_bfs[ i ]['_memo']['picked'] = True

        self.leaf_set = set([ n for n in self.node_set if n.features['_memo']['picked'] and not n.features['_memo']['_children']])
        node_set = kNodes()
        ndel = 0
        for n in self.node_set :
            if n.features['_memo']['picked']:
                n.children = n.features['_memo']['_children']
                n.is_a = n.features['_memo']['_is_a']
                node_set.add(n)

            else:
                ndel +=1
            #del(n.features['_memo'])
        self.node_set = node_set
        print (f"Minimal tree dimensions : {len(self.node_set)}{len(self.leaf_set)}")
        return ndel

    def preload_node(self, ID, name, ora_data, child_node_ids, uniprot_ids, is_leaf=False):
        # create or fetch from pool
        maybe_new_node = OraNode(ID, name, uniprot_ids)
        new_node       = self.node_set.add(maybe_new_node)
        if is_leaf:
            self.leaf_set.add(new_node)
        new_node.children = child_node_ids
       
        if not 'fisher' in new_node.features:
            new_node.features['fisher'] = ora_data
        return new_node
    
    def wire(self):
        index = { n.ID : n for n in self.node_set }
        # suppose we only have children ref
        # ie node.is_a = []

        for node in self.node_set:          
            _ = [ index[n_id] for n_id in node.children ]
            node.children = _
            for n in node.children:
                n.update('is_a', node)
                
            #if not len(n.is_a) == len(set(n.is_a)):
            #    raise TypeError(f"RDR {n.name} {[_.name for _ in n.is_a]}")
        for n in self.node_set:
            if not n.is_a:
                self.root.children = [n]
                return
        raise OraTreeWireError("not primary node found")
            #node.is_a     = [ index[n_id] for n_id in node.is_a ]
            #node.children = [ index[n_id] for n_id in node.children ]
            #print(f"Wire {node.name}:\n\t{[ _.ID for _ in node.is_a]}\n\t{[ _.ID for _ in node.children]}")
    @property
    def dimensions(self):
        return ( len(self.node_set), len(self.leaf_set) )

    def dictify(self,ns_prefix): 
        assert(len (self.root.children) == 1)

        if ns_prefix:
            return self.root.children[0].dictify()

        return [ n.dictify() for n in self.root.children[0].children ]
      
    def jsonify(self): 
        return json.dumps(self.dictify(True))
