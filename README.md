# Manipulating Ontologies

## DAG extraction

* Proteome majoritaire
** Proteome majoritaire +  Elements de proteomes minoritaires
*** Union des protéomes


## Micro-services

### CLI
#### ORA analysis
```bash
python -m unigo pwas cli --goport=5555 --exp-prot=data/lst/UP0000001418_random_exp.lst --de-prot=data/lst/UP0000001418_random_delta.lst --vectorized --taxid=710127
```

#### unigo local_test

#### unigo gostore
Microservice providing following ressources

* `/`: index of currently loaded unigo trees
    
* `/taxids`: list of taxids
    
* `/unigo/<taxid>`: fetch taxid unigo as json

#### unigo pwas test

```sh
python -m unigo pwas test --prot=data/proteomes/uniprot-proteome_UP000000625.xml
```

#### unigo pwas api

##### Vectorized
```json
{'children': [
        {'best': 0.41097679695408007,
   'children': [
                {'bkgFreq': 0.0027013752455795677,
     'maxMemberCount': 22,
     'name': 'leucine biosynthetic process',
     'pvalue': 0.41097679695408007,
     'table': [
                        [
                            5,
                            1551
                        ],
                        [
                            28,
                            10665
                        ]
                    ],
     'uniprotID': ['P0A6A6', 'P0AB80', 'P30125', 'P09151', 'P30126'
                    ]
                }
            ],
   'name': 'leucine biosynthetic process_group'
        },
        {'best': 0.41917231526072896,
   'children': [
                {'bkgFreq': 0.00036836935166994106,
     'maxMemberCount': 3,
     'name': 'response to lithium ion',
     'pvalue': 0.41917231526072896,
     'table': [
                        [
                            1,
                            1555
                        ],
                        [
                            3,
                            10694
                        ]
                    ],
     'uniprotID': ['P0AFA7'
                    ]
                },
                {'bkgFreq': 0.00036836935166994106,
     'maxMemberCount': 3,
     'name': 'intracellular pH reduction',
     'pvalue': 0.41917231526072896,
     'table': [
                        [
                            1,
                            1555
                        ],
                        [
                            3,
                            10694
                        ]
                    ],
     'uniprotID': ['P0AFA7'
                    ]
                }
            ],
   'name': 'response to lithium ion_group'
        }
    ],
 'name': 'root'
}
```

#### unigo pwas cli

## Ontologie

Les termes GO sont organisés en une structure hiérarchique ([DAG](https://en.wikipedia.org/wiki/Directed_acyclic_graph)). Dans ce type de graph les liens sont orientés, mais les noeuds fils peuvent avoir plusieurs parents.

Prenez un moment pour vous familliariser avec [sa structure](https://www.ebi.ac.uk/QuickGO/).

Les termes de cette ontologie sont organisés en trois arbres indépendants:
* **biological process**
* **molecular function**
* **cellular component**

On appellera **namespace** chacun de ces trois arbres. On travaillera séparement sur chaque namespace. Nous vous suggerons de commmencer à travailler sur le namespace **biological process**.

### Topologie des arbres de termes GO
Les noeuds sont des termes GO, auxquels nous allons associer des pValue d'enrichissement.

Les feuilles sont les protéines directement annotées par un terme GO.

###### Attention, règle du vrai chemin

*Si une protéine est porteuse d'un terme GO, alors est elle aussi porteuse de tous les parents de ce terme.*

### Comment lire et manipuler une arbre de termes

######  Après le chargement initial de l'ontologie générique  
`go.setOntology(dataDir + "/go.owl")`

##### On créé  l'abre des termes GO membres d'un namespace et qui annotent  une liste de protéines
```python
goTreeObj = go.createGoTree(ns="biological process", proteinList=xpProtList, uniprotCollection=uniprotCollection)
```
Où
* `ns`, est le namespace GO étudié
* `proteinList`, une liste d'identifiants Uniprot
* `uniprotCollection`, une collection d'objet Uniprot



##### API de l'arbre d'annotations


###### Taille de l'abre
Affiche un tuple (feuilles, noeuds, protéines)
```python
goTreeObj.dimensions
```
#### Extraire un sous arbre
```python
sousArbre = goTreeObj.newRoot("transmembrane transport")
```

#### Enumerer les protéines porteuses de termes GO

##### Dans tout l'arbre
```python
proteinList = goTreeObj.getMembers()
```
##### A partir d'un noeud
```python
proteinList = goTreeObj.getMembersByName("transmembrane transport")
proteinList = goTreeObj.getMembersByID("GO:0055085")
```

#### Acceder directement à un noeud
```python
goTerm = goTreeObj.getByName("transmembrane transport")
goTerm = goTreeObj.getByID("GO:0055085")
```

#### Parcourir tous les noeuds d'un arbre
```python
for goTerm in goTreeObj.walk():
    print(goTerm)
```

#### Elager un arbre
On définit une fonction test (retournant vrai ou faux) à appliquer à chaque noeud de l'arbre de départ.
La fonction test recevra en argument le terme GO/objet noeud à analyser.
Si la fonction retourne faux pour un noeud, ce noeud et tous ses descendants seront éliminés dans l'arbre d'arrivée. 

**Exemple**: retirer tous les terme GO n'ayant aucune des protéines `'P75936', 'P76231', 'P0A8S9'` parmi leurs protéines annotées.
```python
def predicat(goTerm):
    return set(goTerm.getMembers()) & set(['P75936', 'P76231', 'P0A8S9'])

goTreeObj_avec_P75936_P76231_P0A8S9 = goTreeObj.drop(predicat)
```

### Chargements des données uniprot et GO génériques
Avec exemple de création d'un arbre de termes GO du namespace "**biological process**".

Afin de vous familliariser, essayer de réaliser les opérations ci-dessus sur l'arbre `goTreeObj`






## Analyse de l'enrichissement en terme GO parmi les protéines surabondantes - suite

### Analyse de la surreprésentation
On pourrait appliquer la fonction `righEnd_pValue` manuellement à chaque terme de l'arbre.
Mais, pour nous faciliter le travail, la fonction **computeORA** du package *stat_utils.py* permet d'appliquer l'analyse de surreprésentation récursivement, en profondeur, à partir d'un terme GO racine.



#### Estimation récursive des enrichissements

* Choix du terme GO parent à partir duquel l'analyse ORA sera recursivement appliquée
```python
pathWayRoot = xpGoTree.getByName("transmembrane transport")
```

* Définition du terme GO regroupant **tout le protéome**
```python
pathWayBKG = fullEcoliGoTree.getByName("biological process")
```

* Calcul de l'enrichissement en termes GO successifs parmi les protéines surabondantes (ici, *saList*)
```python
oraScores = computeORA(pathWayRoot, saList, pathWayBKG)
```

* Des arbres peuvent également être passés.
```python
tm_transport_GoTree = xpGoTree.newRoot("transmembrane transport")
oraScores = computeORA(tm_transport_GoTree, saList, fullEcoliGoTree)
```

#### Accès aux score d'enrichissements
Vous avez deux moyens d'accéder aux pvalue des différents termes GO

##### La variable oraScores

Elle contient une liste de tuples (pvalue, termeGO), ainsi `print(oraScores[0])` affichera

```python
(0.3826117782409264,
{
    "ID": "GO: 0055085",
    "name": "transmembrane transport",
    "eTag": [
        "P37624",
        "P0AGH1",
        "P75797",
        "P31550",
        "P77348",
        "P06149",
        "P0AF98",
        "P76397",
        "P0AFH2",
        "P0AFH6",
        "Q47622",
        "P60778",
        "P76185",
        "P23843",
        "P0A9V1",
        "Q46863",
        "P77338"
    ],
    "leafCount": 0,
    "features": {
        "Fisher": 0.3826117782409264,
        "Hpg": 0.38143716334525035
    },
    "oNode": obo.GO_0055085,
    "isDAGelem": True,
    "children": [
        "ion transmembrane transport",
        "organic acid transmembrane transport",
        "import across plasma membrane",
        "protein transmembrane transport",
        "carbohydrate transmembrane transport",
        "drug transmembrane transport",
        "purine-containing compound transmembrane transport",
        "nucleoside transmembrane transport",
        "export across plasma membrane",
        "oligopeptide transmembrane transport",
        "regulation of transmembrane transport",
        "transmembrane transporter activity",
        "negative regulation of transmembrane transport",
        "thiamine transmembrane transport",
        "pyrimidine-containing compound transmembrane transport",
        "positive regulation of transmembrane transport",
        "polyamine transmembrane transport"
    ]
})
```
##### L'arbre GO des protéines experimentales 
Chaque noeud stocke sa pvalue, ainsi le parcourt de **l'arbre des protéines expérimentales** suivant 
```python
for n in pathWay.walk():
    print( n.pvalue, n.name, len(n.getMembers(nr=True)) )
```
affichera
```
0.3826117782409264 transmembrane transport 94
0.026319852884717568 ion transmembrane transport 44
None cation transmembrane transport 28
None inorganic cation transmembrane transport 18
None proton transmembrane transport 14
None electron transport coupled proton transport 4
None ATP synthesis coupled proton transport 7
None plasma membrane ATP synthesis coupled proton transport 2
None copper ion transmembrane transport 1
None copper ion export 1
```


#### Mise en forme des résultats des enrichissements
Il vous est demander de créer la fonction `printRankings` qui affichera les pathways par ordre croissant de p-value.
Elle affichera pour chaque terme:
* Le nom du pathway
* La pvalue
* Le nombre de protéines surabondantes membres de ce pathway


#### Représentation des résultats
Afin de faciliter le travail des séances suivantes, vous devez produire une représentation des résultats que vous venez de construire dans un dictionnaire de la forme.
```json
    { 
        "plasma membrane" : {
            "name" : "plasma membrane",
            "pvalue" : 0.9999980015962491,
            "proteineTotal" : [uniprotID, ..],
            "proteineSA" : [uniprotID, ..]
    
           },
       ....
     }
```
Puis vous l'écrirez dans un fichier au [format JSON](https://en.wikipedia.org/wiki/JSON).

#### Représentation de l'objet arbre
On pourrait souhaiter travailler sur la représentation graphique de l'analyse d'enrichissement en s'appuyant sur l'ontologie. On a donc besoin d'ecrire l'état de nos objets **GoTree** pour pouvoir les réulitlier ulterieurement (lors du TP5). La librarie standard fournit le module [pickle](https://docs.python.org/3/library/pickle.html). Cette librarie a deux avantages comparés à la sérialisation en JSON:
* La description de l'objet est automatique, vous n'avez pas besoin d'explicitement spécifier les attributs à sauvgegarder.
* Les méthodes sont aussi sauvegardées, à la désérialieation la variable produite est donc un objet de "plein-droit"

Le format Pickle a cependant quelques limites, essayer d'ecrire un fichier au format pickle représentant l'objet **goTree** et prêtez attention au message s'affichant.

python -m unigo pwas --p-port 1234
