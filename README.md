# Guide Utilisateur

## Sommaire :  <a id="sommaire"></a>

[- Installation et démarrage de Flask](#flask)


## Composition du site
Le site (http://sentinellesbioinfo.pythonanywhere.com) est composé de plusieurs pages dont chacune d'entre-elles a une fonction particulière.

![Zozor](https://zupimages.net/up/20/50/ka5n.png)

### La page d'accueil
Sur cette page, on retrouve les informations sur le thème du site dans un paragraphe "A propos" ainsi qu'un résumé des outils mis à disposition pour l'utilisateur.

![Zozor](https://zupimages.net/up/20/50/jvdk.png)


### Infos générales
Cette page nous donne des renseignements sur l'espèce qui est représenté sur le site : la vipère péliade.
On y trouve une rapide présentation de l'espèce avec :

* ces caractères morphologiques et biologiques,
* des indications sur son comportement et son alimentation,
* sa phylogénie avec sa taxonomie et l'arbre phylogénétique de son espèce.

![Zozor](https://zupimages.net/up/20/50/7uxv.png)

### Séquences

1. l'utilisateur doit choisir les séquences qu'il veut étudier en cochant les cases correspondantes
2. récupération des séquences choisies dans la base de données en ligne (NCBI)

![Zozor](https://zupimages.net/up/20/50/t1qr.png)


### Alignements
Sur cette page, l'utilisateur peut choisir :

1. la méthode d'alignement : Clustal ou Muscle
2. la méthode pour construire l'arbre phylogénétique : Neighbor Joining ou Maximum Likelihood

![Zozor](https://zupimages.net/up/20/50/juk3.png)

3. pour générer l'arbre, l'utilisateur doit cliquer sur "Faire l'arbre"

Mais il peut aussi télécharger un fichier fasta qui contient les séquences choisies précédemment.
![Zozor](https://zupimages.net/up/20/50/awud.png)


### Phylogénie
Cette page donne le résultat de l'alignement sous forme d'arbre.

![Zozor](https://zupimages.net/up/20/50/cb9i.png)

L'utilisateur peut :

* Télécharger le fichier d'alignement précédemment créé

![Zozor](https://zupimages.net/up/20/50/gx41.png)

* Enregistrer l'image de l'arbre en faisant un clic-droit sur l'image

![Zozor](https://zupimages.net/up/20/50/mayf.png)

* Récupérer l'arbre généré en format newick

![Zozor](https://zupimages.net/up/20/50/aos8.png)

-------------------------------------------

[Retourner au sommaire](#sommaire)

-------------------------------------------

# Installation, configuration et démarrage du serveur Flask <a id="flask"></a>

1. Installer Python 3

https://www.python.org/downloads/

2. Cloner le repository GitHub ou télécharger le `.zip`

https://github.com/les-sentinelles-de-la-bioinfo/ProjetPhylogenie

3. Créer l'environnement virtuel dans le répertoire contenant `flask_app.py`

`python3 -m venv env`

4. Activer l'environnement virtuel

`source env/bin/activate`

5. Installation de Flask et des modules/packages
```
pip3 install Flask
pip3 install matplotlib
pip3 install biopython
```

6. Démarrage de Flask

```
export FLASK_APP=flask_app.py
export FLASK_ENV=env
flask run
```

7. L'interface utilisateur est maintenant disponible en local à l'adresse suivante : http://127.0.0.1:5000

----------------

[Retourner au sommaire](#sommaire)
