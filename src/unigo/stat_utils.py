"""
Basic statistical functions to compute pathway enrichment

"""

from scipy.stats import hypergeom
from scipy.stats import fisher_exact
import numpy as np
from scipy.cluster.hierarchy import dendrogram, ward, fcluster

def computeSelfORA(node, proteinList):
    """computeSelfORA
    Compute the probability of observing at least the k members of the supplied GOnodes/pathway
    present in the supplied protein list
    
    Parameters
    ----------
    node : element of the ontology tree.
    proteinList: list or uniprot identifiers
    
    Remarks
    ----------
    Background proteins are all uniprot ID found in the supplied annotation subtree
    """

    ORA = []
    
    # universe is all uniprotID found in the annotation tree
    universe = set(node.getMembers()) 
    N = len(universe)
    # nSet is the observed set   
    nSet = set(proteinList)
    n = len(nSet)
    for cPath in node.walk():       
        Kstates = set(cPath.getMembers())
        K = len( Kstates )
        print(f"{cPath.name} has {K} members")
        
        k_obs = Kstates & nSet
        k = len(k_obs)
        
        p = righEnd_pValue(N, n, K, k)
        
        ORA.append( (p, cPath) )
        print(f"{cPath.name} [{K} -> {k}/{n}] = {p}")
    return ORA
    
def righEnd_pValue(N, n, K, k):
    
#print(f"N={N}, n={n}, K={K}, k={k}")

#The hypergeometric distribution models drawing objects from a bin. 
#N is the total number of objects, K is total number of Type I objects. 
#The random variate represents the number of Type I objects in N 
#drawn without replacement from the total population.

# Right-end tail of the CDF is P(X>=k)
    p_x = hypergeom(N, K, n).cdf([k - 1])
    return 1.00 - p_x[0]

def computeORA(node, proteinList, nodeBKG, verbose=False):
    Fisher, CDF = computeORA_BKG(node, proteinList, nodeBKG, verbose=verbose)
    return Fisher

def computeORA_BKG(node, proteinList, nodeBKG, verbose=False): # IDEM, mais avec un autre arbre de reference
    
    if verbose:
        print("Computing Over Representation w/ a background tree")

    ORA_Fisher = []
    ORA_CDF = []

    universe = set(nodeBKG.getMembers())
    o = len(universe)

    # nSet is the observed set   
    nSet = set(proteinList)
    n = len(nSet)
    
    pathwayPotential = 0
    pathwayReal = 0
    
    import time

    start = time.time()
   # for cPath in node.walk(mustContain=proteinList):
    for cPath in node.walk():
        pathwayPotential += 1
        #verbose = cPath.name == 'enzyme binding'
        if verbose:
            print(cPath.name)

        # Table de contingence
        #
        #        | Pa  | non_PA |
        # -----------------------
        #    SA  |     |        |
        #  nonSA |     |        |

        
        # l'intersection entre les protéines porteuse de l'annotation courante et 
        # la liste des protéines sur-abondante
        # => nombre de succès observés, "k"
        Kstates = set(cPath.getMembers())
        k_obs = Kstates & nSet
        if not k_obs:
            if verbose:
                print("k_obs == 0")
            continue
        k = len(k_obs)
        pathwayReal += 1
        
        # Pour estimer le nombre de protéines non surAbondantes appartenant au pathway ou non
        # Nous utilisons la proporition de protéines du pathway ou non dans le protéome entier
        bkgPath = nodeBKG.getByName(cPath.name)
        bkgPathFreq = len( set(bkgPath.getMembers()) ) / len(universe)  # Fraction du protéomes appartenant à ce Pathway
        nSA_Pa = int ( (o - k) * bkgPathFreq )
        nSA_nPa = int( (o - k) - nSA_Pa )
        
        TC = [
            [ k ,  len(proteinList) - k],
            [ nSA_Pa ,  nSA_nPa]
        ]
        
        oddsratio, pValue = fisher_exact(TC, alternative="greater")
        p = righEnd_pValue(o, n, len( set(bkgPath.getMembers()) ), k)

        if verbose:
            print(f"{cPath.name} {TC} p={pValue} // pL={p}")

        
        ORA_Fisher.append( (pValue, cPath) ) 
        ORA_CDF.append( ( p, cPath) ) 
        
        cPath.set(Fisher=pValue, Hpg=p)

    end = time.time()    
    print(f"Evaluated {pathwayReal} / {pathwayPotential} Pathways, based on {n} proteins in {end - start} sc")
    return ORA_Fisher, ORA_CDF

def applyOraToVector(vectorizedProteomeTree, experimentalProteinID, deltaProteinID, threshold=0.005):
    """ compute Fischer exact test on a list of datastructure corresponding 
        to a vectorized full proteome GO tree.

        Parameters
        ----------
        vectorizedProteomeTree : Univgo.vectorize() return datastructure
        experimentalProteinID  : List of uniprot identifiers of observed proteins
        deltaProteinID         : List of uniprot identifiers of over/under expressed proteins
    """

    def ora(universe, pathway, observedProteinList, deltaProteinList):
        o       = len(universe)
        nSet    = set( observedProteinList )
        n       = len(nSet)
        Kstates = set( pathway["elements"] )
        k_obs   = Kstates & set(deltaProteinList)
        k = len(k_obs)
        
        nSA_Pa  = int ( (o - k) *  pathway["freq"] )
        nSA_nPa = int( (o - k) - nSA_Pa )
            
        TC = [
            [ k ,  n - k],
            [ nSA_Pa ,  nSA_nPa]
        ]

        oddsratio, pValue = fisher_exact(TC, alternative="greater")
        return {
            "name"       : pathway["name"],
            "pvalue"     : pValue,
            "K_states"   : pathway["elements"],
            "k_success"  : list(k_obs),
            "table"      : TC,
            "bkgFreq"    : pathway["freq"]
        }

    d = vectorizedProteomeTree
    registry = d["registry"]
    # delta _include_in experimental _include_in wholeProteome
    assert( not set(deltaProteinID)        - set(experimentalProteinID) )
    assert( not set(experimentalProteinID) - set(registry)              )
    # Convert two uniprotID list to integers
    expUniprotIndex   = [ d["registry"].index(_) for _ in experimentalProteinID ]
    deltaUniprotIndex = [ d["registry"].index(_) for _ in deltaProteinID        ]
    
    res = { goID: ora(d['registry'], ptw, expUniprotIndex, deltaUniprotIndex) for goID, ptw in d['terms'].items() }

    return { k:v for k,v in res.items() if v["pvalue"] < threshold }




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
    print(V)

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
