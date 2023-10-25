from .tree import GO_tree
from uniprot_redis.store.schemas import UniprotDatum
from typing import Literal
"""
We wanna reduce sub tree crawling operations
To apply the "one protein to a node -> one protein to all parent nodes" rule
We start from leaves of annotation tree and ascend

"""

type SortCrit = Literal['pvalue', 'count', 'bkfq']

class ORA_tree(nx.DiGraph):
    def __init__(self, src_tree:GO_tree, xp_uniprots:Iterator[UniprotDatum]): # we delegated to external obolete uniprot id forwarding
        super().__init__()
        assert(src_tree.ora_rdy)
        def go_up(base_tree, node_id, )
        for go_id in src_tree.leave_ids:
          
            n_dict              = src_tree.get_go_node(go_id)
            uniprot_data_buffer = n_dict['proteins'] if 'proteins' in n_dict else set()

            self.add_node(go_id, **n_dict)


    def sorted_by(self, crit:SortCrit='pvalue'):
        pass



