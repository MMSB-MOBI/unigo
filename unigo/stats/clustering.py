
def kappaClustering(registry, applyOraToVectorResults, fuseThresh=0.2):
    """ Cluster Pathways evaluated by applyOraToVector scoring using Cohen's kappa coefficient
    """
    omega = set( [ k for k in range(len(registry)) ] )
    pathwayID = list( applyOraToVectorResults.keys() )

    for x in range(0, len(pathwayID)):
        _ = applyOraToVectorResults[ pathwayID[x] ]
    #    print(x," | ", _["name"], _["K_states"])



    #pathwayID = pathwayID + pathwayID
    pdist = []  
    for i in range(0, len(pathwayID) - 1):
        iID, iTerm   = (pathwayID[i],  applyOraToVectorResults[ pathwayID[i] ])
        iName        = iTerm["name"]

        in_iTerm     = set(iTerm["K_states"])
        not_in_iTerm = omega - in_iTerm
        for j in range(i + 1, len(pathwayID)):
            jID, jTerm = (pathwayID[j],  applyOraToVectorResults[ pathwayID[j] ])
            jName        = jTerm["name"]

            in_jTerm     = set(jTerm["K_states"])
            not_in_jTerm = omega - in_jTerm
    #        print(f"\n\nGO Term [{i}]{iName} [{j}]{jName}")
    #        print(in_iTerm, in_jTerm)#, not_in_iTerm)
            k = kappa( in_iTerm     & in_jTerm    ,\
                       in_iTerm     & not_in_jTerm,\
                       not_in_iTerm & in_jTerm    ,\
                       not_in_iTerm & not_in_jTerm ,\
                )
            
    #        print(f"pdist = {k}")
            pdist.append( 1 - k)

    Z = ward(np.array(pdist))
    #print(Z)
    _Z = fcluster(Z, t=fuseThresh, criterion='distance')
    #print(_Z)

    V = flattenToD3hierarchy(_Z.tolist(), registry, applyOraToVectorResults, pathwayID)
    #print(V)

    return V

def kappa(a, b, c, d):
    """           GO term 2
               | yes |  no |        
    -------------------------------   
    GO   | yes |  a  |  b  | 
   term1 |  no |  c  |  d  |


   kapa(GO_1, GO_2) = 1 - (1 - po) / (1 - pe)

   po = (a + d) / (a + b + c + d) 
   marginal_a = ( (a + b) * ( a + c )) / (a + b + c + d)
   marginal_b = ( (c + d) * ( b + d )) / (a + b + c + d)
   pe = (marginal_a + marginal_b) / (a + b + c + d)

"""
    a = float(len(a))
    b = float(len(b))
    c = float(len(c))
    d = float(len(d))

    po = (a + d) / (a + b + c + d) 
    marginal_a = ( (a + b) * ( a + c )) / (a + b + c + d)
    marginal_b = ( (c + d) * ( b + d )) / (a + b + c + d)
    pe = (marginal_a + marginal_b) / (a + b + c + d)

    #print (f" {a} | {b}\n {c} | {d}")
    return 1 - (1 - po) / (1 - pe)


def flattenToD3hierarchy(_Z, registry, applyOraToVectorResults, pathwayID):
    clusterElement = { }
    for pathwayIndex, clusterNum in enumerate(_Z):
        currPathwayID    = pathwayID[pathwayIndex]
        currPathway      = applyOraToVectorResults[currPathwayID]
        if not clusterNum in clusterElement:
            clusterElement[clusterNum] = {
                "name": None,
                "children":[],
                "best" : 1.1                
            }
        currCluster = clusterElement[clusterNum]
        proxy = { k : v for k,v in currPathway.items() if not (k == "K_states" or  k == "k_success") }
        currCluster["children"].append(proxy)
        proxy["uniprotID"] = [ registry[int(_uid)] for _uid in currPathway["k_success" ] ]
        proxy["maxMemberCount"] = len(currPathway["K_states"])
        if currPathway["pvalue"] < currCluster["best"]:
            currCluster["best"] = currPathway["pvalue"]
            currCluster["name"] = currPathway["name"] + "_group"

    return {
        "name" : "root",
        "children" : [ clust for _, clust in clusterElement.items() ]
    }
