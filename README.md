# Guide Utilisateur

Sommaire :


[- Guide installation et de configuration de WAMP sur Windows](#guidewamp)

[- Guide installation et de configuration de LAMP sur Linux](#guidelinux)


## Composition du site
Le site est composé de plusieurs pages dont chacune d'entre-elles a une fonction particulière.

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

# Site Web sur serveur local WAMP  <a id="guidewamp"></a>
![Zozor](https://blog.nicolashachet.com/wp-content/uploads/2011/05/wamp.png)


## Téléchargement et installation du logiciel Wampserver
https://www.wampserver.com/en/download-wampserver-64bits/


## Insatllation package mod_wsgi

Dans un terminal, tapez les commandes suivantes:
1. pip install mod_wsgi
2. mod_wsgi-express install-module
3. Allez dans votre dossier wamp64\bin\apache\apache2.4.46\conf\httpd.conf
4. Copiez les 3 lignes obtenues en point 2 à la suite des "LoadModule" (l'endroit n'a pas d'importance)

Vous devriez avoir ce genre de lignes

```
# WSGI module
LoadFile "C:/Python37/python37.dll"
LoadModule wsgi_module "C:/Python37/lib/site-packages/mod_wsgi/server/mod_wsgi.cp37-win_amd64.pyd"
WSGIPythonHome "C:/Python37"
WSGIApplicationGroup %{GLOBAL}
````

5. Sauvegardez et fermez le fichier

## Ajout du script .wsgi
Le plus pratique est de mettre ce script dans un dossier dédié. Ici, nous l'appelons, wsgi_scripts.

1. Créez un fichier phylogenie.wsgi et ajoutez ces lignes de code

```
import sys 

#Expand Python classes path with your app's path
sys.path.insert(0, 'C:/wamp64/www/Phylogenie') 

from flask_app import app as application

```

2. Sauvegardez et quittez

## Configuration du virtual host de Wampserver

1. Créez le virtual host depuis la page d'accueil de Wamp (localhost)

Ajoutez un nouvel virtual host (1)
  
![Zozor](https://zupimages.net/up/20/50/fzu4.png)

    * Entrez le nom du fichier (1)
    * Entrez le chemin complet du fichier (2)
    * Lancez la procédure de création (3)

![Zozor](https://zupimages.net/up/20/50/nzdp.png)


2. Allez dans le fichier wamp64\bin\apache\apache2.4.46\conf\extra\httpd-vhosts.conf
3. Modifiez le code de votre virtual host comme ceci (les chemins sont à modifier en fonction de votre configuration) :

```
#
<VirtualHost *:80>
  ServerName phylogenie
  ServerAlias PY
  DocumentRoot "c:/wamp64/www/phylogenie"
  WSGIScriptAlias / "c:/wamp64/www/phylogenie/wsgi_scripts/phylogenie.wsgi"
  <Directory "c:/wamp64/www/phylogenie/">
    Options +Indexes +Includes +FollowSymLinks +MultiViews
    AllowOverride All
    Require local
  </Directory>
</VirtualHost>

```

4. Sauvegardez et fermez le fichier

## Démarrage de Wampserver 

1. Double cliquer sur l'icone du logiciel
2. Quand le serveur est prêt, il apparait en vert dans le sous-menu de la barre des tâches

![Zozor](https://zupimages.net/up/20/50/0jwl.bmp)

3. Apuyez dessus avec le clic gauche et choisiez "Vos VirtualsHosts"

![Zozor](https://zupimages.net/up/20/50/zygb.bmp)

4. Cliquez sur celui que vous avez créé et qui apparait dans la liste

![Zozor](https://zupimages.net/up/20/50/ehzt.bmp)

-----------------------------
# Configuration utilisée

* Windows 10 64-bits 20H2
* Python 3.7
* biopython 1.78
* numpy 1.19.3
* mod_wsgi 4.7.1
* matplotlib 3.3.3
* Flask 1.1.2 

----------------------------

----------------------------

# Site Web sur serveur Apache Linux  <a id="guidelinux"></a>

Les fichiers sont à placer dans /var/www/html/

-------------------

## Installation des librairies                 


Biopython : `pip install biopython`

Mettre à jour Biopython : `pip install biopython --upgrade`

Flask : `pip insatll Flask`

matplotlib : `pip install matplotlib`

------------------

1. Installation de apache2

`sudo apt install apache2 php libapache2-mod-php mysql-server php-mysql`

2. Installation du mod_wsgi

`pip install mod_wsgi-standalone`

3. Redémarrage de Apache

`sudo systemctl restart apache2`

4. Modification du virtualhost

s`udo nano /etc/apache2/conf-available/mod-wsgi.conf`

5. Ajouter la ligne suivante

`WSGIScriptAlias /Phylogenie /var/www/html/wsgi_scripts/phylogenie.wsgi`

6. Appliquer la nouvelle configuration

`sudo a2enconf mod-wsgi`

7. Redémarrer Apache

`sudo systemctl restart apache2`

8. Créez un dossier wsgi_scripts et un fichier phylogenie.wsgi et ajoutez ces lignes de code

```
import sys 

#Expand Python classes path with your app's path
sys.path.insert(0, 'var/www/html/Phylogenie') 

from flask_app import app as application

```

9. Configuration du virtualhost de Flask

`sudo nano /etc/apache2/sites-available/FlaskApp.conf`

10. Taper les lignes suivantes

```
<VirtualHost *:80>
	ServerName localhost
	WSGIScriptAlias / /var/www/html/Phylogenie/wsgi_scripts/phylogenie.wsgi
	<Directory /var/www/html/Phylogenie/>
		Order allow,deny
		Allow from all
	</Directory>
	Alias /static /var/www/html/Phylogenie/static
	<Directory /var/www/html/Phylogenie/static/>
		Order allow,deny
		Allow from all
	</Directory>
	ErrorLog ${APACHE_LOG_DIR}/error.log
	LogLevel warn
	CustomLog ${APACHE_LOG_DIR}/access.log combined
</VirtualHost>
```
Sauver et quitter

11. Appliquer la nouvelle configuration

`sudo a2ensite FlaskApp`

12. Redémarer Apache

`sudo service apache2 restart`

Le site est maintenant accessible à l'adresse http://localhost/Phylogenie 
