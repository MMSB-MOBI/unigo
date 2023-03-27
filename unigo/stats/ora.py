"""
Basic statistical functions to compute pathway enrichment

"""
from scipy.stats import hypergeom
from scipy.stats import fisher_exact
from ..tree.ora_results import OraTree

def righEnd_pValue(N, n, K, k):
    
#print(f"N={N}, n={n}, K={K}, k={k}")

#The hypergeometric distribution models drawing objects from a bin. 
#N is the total number of objects, K is total number of Type I objects. 
#The random variate represents the number of Type I objects in N 
#drawn without replacement from the total population.

# Right-end tail of the CDF is P(X>=k)
    p_x = hypergeom(N, K, n).cdf([k - 1])
    return 1.00 - p_x[0]



def node_ora(node, SA, tot_pop_deep, tot_pop_annot, verbose):
    # contingency table
    #
    #        | Pa  | non_PA |
    # -----------------------
    #    SA  | (A) |   (B)  |
    #  nonSA | (C) |   (D)  |

    card_tot_pop_deep = len(tot_pop_deep)
    card_tot_pop_annot = len(tot_pop_annot)
    # Pathway protein set
    Pa_deep       = set(node.background_members)
    Pa_annot      = set(node.getMembers())
    Pa_deep_card  = len(Pa_deep)
    Pa_annot_card = len(Pa_annot)
    
#    print(f"Pa_deep_card:{Pa_deep_card} , tot_pop:{card_tot_pop_deep}")

 #   print(f"SA & Omega_deep {len(SA & tot_pop_deep)}")
 #   print(f"Omega_deep - SA {len(tot_pop_deep - SA)}")
    #  A) SA & Pa intersection
    #  aka observed succes states "k"
    assert len ( (Pa_deep & SA) - (Pa_annot & SA) ) == 0
    SA_Pa = Pa_annot & SA
    if not SA_Pa:
        if verbose:
            print("k_obs == 0")
        return None
    k = len(SA_Pa)
    
    # B) SA & nPa intersection
    # Any protein overabundant and not member of the pathway
    SA_nPA_deep = SA - Pa_deep
    card_B_deep = len(SA_nPA_deep)
    SA_nPA_annot = SA - Pa_annot
    card_B_annot = len(SA_nPA_annot)
    # C) nSA & Pa
    # Any protein not overabundant and member of the pathway
    nSA_Pa_deep = Pa_deep - SA
    card_C_deep = len(nSA_Pa_deep)
    nSA_Pa_annot = Pa_annot - SA
    card_C_annot = len(nSA_Pa_annot)
    
    # D) nSA & nPA
    # It is the remaining sub population 
    #(D) = card_tot_pot - [Card(A) + Card(B) + Card(C)]
    # DEBUG
    #D_deep = tot_pop_deep -  SA_Pa - SA_nPA_deep - nSA_Pa_deep
    #print(f"Explicit Card_D is {len(D_deep)}")
    #card_D_deep = len(D_deep)
    card_D_deep  = card_tot_pop_deep  - (k + card_B_deep + card_C_deep)
    card_D_annot = card_tot_pop_annot - (k + card_B_annot + card_C_annot)
    TC_deep = [
        [ k ,  card_B_deep],
        [ card_C_deep,  card_D_deep]
    ]
  #  print(TC_deep)
   
    oddsratio, pValue_deep = fisher_exact(TC_deep, alternative="greater")
    TC_annot = [
        [ k,  card_B_annot],
        [ card_C_annot,  card_D_annot]
    ]
   # print(TC_annot)
    oddsratio, pValue_annot = fisher_exact(TC_annot, alternative="greater")
    #p = righEnd_pValue(card_tot_pot, len(SA), len(Pa), k)

    if verbose:
        print(f"{node.name} {TC_deep} p={pValue_deep}")
    return {
            "name"              : node.name,
            "ID"                : node.ID,
            "pvalue_deep"       : pValue_deep,
            "pvalue_annot"      : pValue_annot,
            "K_k_N_n_deep"      : (Pa_deep_card, k,  card_tot_pop_deep, len(SA)),
            "K_k_N_n_annot"     : (Pa_annot_card, k, card_tot_pop_annot, len(SA)),
            "table_deep"        : TC_deep,
            "table_annot"       : TC_annot,
            "xp_hits"           : list(Pa_annot & SA),
            "xp_obs"            : list(Pa_annot),
            "pathway_freq_deep" : node.background_frequency,
            "pathway_freq_obs"  : len(Pa_annot & SA) / card_tot_pop_annot, 
        }

"""
ORA Computation over a masked unigo wrapper
We order by pvalue_dee for now, to discuss
"""
def apply_ora_to_unigo(unigo_obs_masked, sa_protein_list, verbose=False, pvalue_max = 0.1):
    if not unigo_obs_masked.masked:
        raise AttributeError("Applying ORA compuation on a non-masked __aka blueprint__ tree is forbidden")
    # We must remove from the sa_protein_list the protein with no annotation in current tree
    sa_protein_list = list(set(sa_protein_list) & unigo_obs_masked.tree.proteins)
    ora_res_Fisher = _apply_ora_to_tree(unigo_obs_masked.tree, \
                                             sa_protein_list, verbose=verbose)
    
    return [ v for (k,v) in sorted(ora_res_Fisher.items(), key=lambda item: item[1]["pvalue_deep"]) \
            if v["pvalue_deep"] <= pvalue_max ]

"""
Take a unigo, a set of observed and delta
Population Background can be set the blueprint population __aka proteome__ or the 
current set of observed protein
Nota:
Tree nodes are currently modified to carry pvalues for no later purpose yet,
it is unecessaary and kinda hacky wayto prserve GO term organisation
"""
def _apply_ora_to_tree(tree, proteinList, verbose=False): # IDEM, mais avec un autre arbre de reference
    
    if verbose:
        print("Computing Over Representation w/ a background tree")

    node = tree.root
    # The SurAbundant set is invariant (nb of draw), so is the total population (nb_states_total)
    SA = set(proteinList)
    tot_pop_deep = set(node.background_members)
    tot_pop_annot = set(node.getMembers())
    
    pathwayPotential = 0
    pathwayReal = 0
    
    import time

    start = time.time()
    
    ora_res = {}
    for cPath in node.walk():
        pathwayPotential += 1
        ora_node = node_ora(cPath, SA, tot_pop_deep, tot_pop_annot, verbose)
        if ora_node is None:
            continue
        pathwayReal += 1    
        ora_res[cPath.name] = ora_node
        cPath.set(Fisher=ora_node)

    end = time.time()    
    print(f"Evaluated {pathwayReal} / {pathwayPotential} Pathways, based on {len(SA)} proteins in {end - start} sc")
    
    return ora_res  
  

# Reflect masked unigo tree in a ora score tree
# then shove non significant term from bottom
def apply_ora_tree2tree(unigo_obs_masked, unk_uniprot_ids, sa_protein_list, ref_pop:str, verbose=False)->OraTree:

    tree = unigo_obs_masked.tree
    # We merge protein with no annotation with the background population for Fisher table consistency
    tot_pop_deep  = set(tree.root.background_members) | unk_uniprot_ids
    tot_pop_annot = set(tree.root.getMembers()) | unk_uniprot_ids
    # The SurAbundant set is invariant (nb of draw), so is the total population (nb_states_total)
    SA = set(sa_protein_list)
    print(f"apply_ora_tree2tree over {unigo_obs_masked.ns} ontology\n",
          f"number of proteome proteins with annotation        : {len(tree.root.background_members)}\n",
          f"number of experimental proteins with annotation    : {len(tree.root.getMembers())}\n",
          f"number of experimental proteins without annotation : {len(unk_uniprot_ids)}\n",
          f"numer of overabundant proteins                     : {len(SA)}\n",
          f"numer of overabundant proteins  without annotation : {len(SA & unk_uniprot_ids)}\n",
          f"Background pop. size (proteome and exp. observed)  : {len(tot_pop_deep)},{len(tot_pop_annot)}\n"
        )

    ora_tree = OraTree()
    prev_ora_node = ora_tree.root
    for curr_pathway in tree.walk_no_root():
        if not len([ p.ID for p in curr_pathway.children ]) == len(set([ p.ID for p in curr_pathway.children ])):
            print("apply_ora_tree2tree::WARNING: masked unigo pathway RR###", [ p.ID for p in curr_pathway.children ])
        ora_data = node_ora(curr_pathway, SA, tot_pop_deep, tot_pop_annot, verbose)
        curr_ora_node = ora_tree.preload_node(  curr_pathway.ID, curr_pathway.name, 
                                            ora_data, [ p.ID for p in curr_pathway.children ], \
                                            curr_pathway.getMembers(nr=True),
                                            len(curr_pathway.children) == 0)
        prev_ora_node = curr_ora_node
    
    ora_tree.wire()
    return ora_tree

    
def applyOraToVector(vectorizedProteomeTree, experimentalProteinID, deltaProteinID, threshold=0.05, translateID=False):
    """ compute Fischer exact test on a list of datastructure corresponding 
        to a vectorized full proteome GO tree.

        Parameters
        ----------
        vectorizedProteomeTree : Univgo.vectorize() return datastructure
        experimentalProteinID  : List of uniprot identifiers of observed proteins
        deltaProteinID         : List of uniprot identifiers of over/under expressed proteins
    """
    d = vectorizedProteomeTree
    registry = d["registry"]

    def ora_obs(universe, pathway, observedProteinList, deltaProteinList):
        # Table de contingence
        #
        #        | Pa                        | non_PA                      |
        # -----------------------------------------------------------------------------
        #    SA  | delta & PA                |  delta - SA_PA              |  delta
        #  nonSA | PA - PA_SA                |  (obs - PA) - SA_nPA        |  obs - delta
        # -----------------------------------------------------------------------------
        #        |  obs & prot_pathway (PA)  |  obs - PA                   |  obs

        #ANNOT CONTINGENCE : 
        #        | Pa                        | non_PA                      |
        # -----------------------------------------------------------------------------
        #    SA  | delta & PA                |  delta - SA_PA              |  delta
        #  nonSA | (annot * freq) - PA_SA    |  (annot - delta) - nSA_PA   |  annot - delta
        # -----------------------------------------------------------------------------
        #        |                           |                             |  annot


        # On implemente deux tables de contingence, une par rapport aux prot observées et une par rapport aux prot annotées (b "devient" d = les protéines annotées)
        #
        #print("xx>",observedProteinList)
        delta = set(deltaProteinList) #surexpressed proteins
        obs = set(observedProteinList) #observed experience proteins
        prot_pathway = set(pathway["elements"]) #all proteins in pathway
        #print("==>", len(obs))
        PA = obs & prot_pathway #observed proteins in pathway
        SA_PA = delta & PA #surrexpressed proteins in pathway
        SA_nPA = delta - SA_PA #surrexpressed proteins not in pathway
        nSA_PA = PA - SA_PA 
        nSA_nPA = (obs - PA) - SA_nPA

        pathway_freq_obs = len(PA) / len(obs) 

        #ANNOT
        annot = len(universe) #WARNING IT'S NOT INDEXED AS OTHER PROTEINS LIST, so works with cardinal
        annot_nSA = annot - len(delta)
        annot_nPA = annot - len(prot_pathway)
        #annot_nSA_PA = len(prot_pathway - SA_PA)
        annot_nSA_PA = annot * pathway["freq"] - len(SA_PA)
        annot_nSA_nPA = annot_nSA - annot_nSA_PA
        if (annot_nSA_nPA) < 0 : 
            annot_nSA_nPA = 0 #TEMPORARY

        TC = [
            [ len(SA_PA), len (SA_nPA)],
            [ len(nSA_PA), len(nSA_nPA)]
        ]

        TC2 = [
            [ len(SA_PA), len (SA_nPA)],
            [ annot_nSA_PA, annot_nSA_nPA]
        ]
    

        oddsratio_obs, pValue_obs = fisher_exact(TC, alternative="greater")
        oddsratio_annot, pValue_annot = fisher_exact(TC2, alternative="greater")
        #print("pvalues", pValue_obs, pValue_annot)
        #print("pathways", pathway["freq"], pathway_freq_obs)
        #print(pathway.keys())
        return {
            "name"       : pathway["name"],
            "pvalue"     : pValue_obs,
            "pvalue_annot" : pValue_annot,
            "K_states"   : pathway["elements"],
            "k_success"  : list(SA_PA),
            "table"      : TC,
            "bkFreq" : pathway["freq"], 
            #"ns" : pathway["ns"], 
            "xp_hits"    : [ registry[iNum] for iNum in SA_PA ] if translateID else None,
            "pathway_freq_annot" : pathway["freq"],
            "pathway_freq_obs" : pathway_freq_obs, 
        }

    assert( not set(deltaProteinID)        - set(experimentalProteinID) )
    annotatedExperimentalProteinID = list( set(experimentalProteinID) & set(registry) )
    annotatedDeltaProteinID        = list( set(deltaProteinID)        & set(registry) )

    if not annotatedExperimentalProteinID:
        msg = "The uniprot ID you supplied could not match "   +\
            "any from the reference proteome/go tree|vector\n" +\
            "Query IDs:" + str(experimentalProteinID) + "\n"   +\
            "Ref IDs:" + str(registry)
        raise ValueError(msg)

    expUniprotIndex   = [ d["registry"].index(_) for _ in annotatedExperimentalProteinID ]
    deltaUniprotIndex = [ d["registry"].index(_) for _ in annotatedDeltaProteinID        ]

    res = { goID: ora_obs(d['registry'], ptw, expUniprotIndex, deltaUniprotIndex) for goID, ptw in d['terms'].items()}
    
    return { k:v for k,v in sorted(res.items(), key=lambda item: item[1]["pvalue"]) if v["pvalue"] < threshold }

