import numpy as np
from scipy.cluster.hierarchy import ward, to_tree

class Tree:
    def __init__(self, distance_matrix, go_terms, scipy_tree):
        self.dist_matrix = distance_matrix
        self.terms = go_terms
        self.root = self.construct_tree(scipy_tree)
        self._name_all_nodes()
    
    def construct_tree(self, scipy_root_node):
        print("# Construct tree")
        root_node = TreeNode(scipy_root_node.id, None)
        self._iterate_to_construct(scipy_root_node, root_node)
        return root_node

    def _iterate_to_construct(self, scipy_node, homemade_node):
        if scipy_node.count == 1: #It's a leaf, give go to it
            _propagate_leaves(homemade_node, homemade_node)
            homemade_node.go = self.terms[scipy_node.id]
            return
   
        new_node_left = TreeNode(scipy_node.left.id, homemade_node)
        new_node_right = TreeNode(scipy_node.right.id, homemade_node)
        homemade_node.children.append(new_node_left)
        homemade_node.children.append(new_node_right)
        self._iterate_to_construct(scipy_node.left, new_node_left)
        self._iterate_to_construct(scipy_node.right, new_node_right)

    def _name_all_nodes(self):
        print("# Name all nodes")
        self._name_node(self.root)

    def _name_node(self, node):
        node.name_it(self.dist_matrix, self.terms)
        if node.children:
            for child in node.children:
                self._name_node(child)

    def serialize_to_d3(self):
        obj = {'name': self.root.name, 'idx': self.root.idx}
        _iterate_serial(self.root, obj)
        return obj

def _iterate_serial(node, serialized):
    if node.children:
        serialized['children'] = []
        for child in node.children:
            new_obj = {'name': child.name, 'idx': child.idx}
            serialized['children'].append(new_obj)
            _iterate_serial(child, new_obj)
    

    
class TreeNode:
    def __init__(self, idx, parent):
        self.idx = idx
        self.parent = parent
        self.children = []
        self.leaves = []
        self.go = None
        self.name = None

    def name_it(self, dist_matrix, terms):
        if not self.leaves:
            raise Exception("You shouldn't have no leaves")
        
        if self.go:
            self.name = self.go['name']
            
        else:
            all_sum_dist = []
            for leaf in self.leaves:
                sum_dist = 0
                for other_leaf in self.leaves:
                    sum_dist += dist_matrix[leaf.idx][other_leaf.idx]
                all_sum_dist.append(sum_dist)
            sum_dist_min = min(all_sum_dist)
            central_leaf = self.leaves[all_sum_dist.index(sum_dist_min)]
            go_term = terms[central_leaf.idx]
            self.name = go_term['name'] + " group"

def create_go_clusters_tree(terms):
    dist_matrix = _distance_matrix(terms)
    Z = ward(dist_matrix)
    tree = to_tree(Z)
    cluster_tree = Tree(dist_matrix, terms, tree)
    return cluster_tree

def _distance_matrix(terms):
    dist_matrix = np.zeros((len(terms),len(terms)))
    for i in range(len(terms)):
        go1 = terms[i]
        for j in range(i+1, len(terms)):
            go2 = terms[j]
            p1 = set(go1["k_success"])
            p2 = set(go2["k_success"])
            dist = len(p1 & p2) / max(len(p1), len(p2))
            dist_matrix[i][j] = dist
            dist_matrix[j][i] = dist
    return dist_matrix

def _browse_tree(scipy_node, homemade_node, terms):
    if scipy_node.count == 1: #It's a leaf, give go to it
        _propagate_leaves(homemade_node, homemade_node)
        homemade_node.go = terms[scipy_node.id]
        return
   
    new_node_left = TreeNode(scipy_node.left.id, homemade_node)
    new_node_right = TreeNode(scipy_node.right.id, homemade_node)
    homemade_node.children.append(new_node_left)
    homemade_node.children.append(new_node_right)
    _browse_tree(scipy_node.left, new_node_left)
    _browse_tree(scipy_node.right, new_node_right)

def _propagate_leaves(node, leaf):
    node.leaves.append(leaf)
    if node.parent:
        _propagate_leaves(node.parent, leaf)

