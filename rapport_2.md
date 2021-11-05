
# Rapport 05/11/2021

## Idées théoriques pour la thèse
Après entretien avec Karen, elle a été convaincue par cette idée :
- Generation procédurale basée sur une liste de missions
	- Implique une future sub task : récuperer une liste de tâches/missions possibles.

L'idée : à partir d'une liste de tâches, un graphe est réalisé pour construire une grande mission, qui sera elle-même traduite en un environnement 3D.
Les "+" : 
- La mission semble pouvoir être lisible avant de voir la simulation
- L'idée peut être relativement simple, intuitive (?)
- Me semble extensible.


Pour cela il faudra utiliser des grammaires (comme L-System par exemple).
- Tache 1 : creer un simple system de grammaire pour les manipuler facilement.
- Tache 2 : Etudier les algos de génération aléatoire (/génétiques/évolutifs) de grammaires.
- Tache 3 : Essayer de voir un premier resultat avec Karen.

On peut utiliser la grammaire spatiale associée à la grammaire de missions pour donner un rendu 3D
- Tache 1 : Essayer une association simple de grammaire spatiale et de mission
- Tache 2 : Voir pour la modélisation 3D avec cette grammaire.
- Tache 3 : Présenter les résultats à Karen.


## Recap du mois :
Premier mois : 
- Reprise de la modelisation 3D opengl.
- Fouille et optimisation de codes simples.
- Création terrain 3D procéduraux par bruit de Perlin -> terrain théoriquement extensible à l'infini.
- Premières tentatives à des algos d'erosion / simulation marins.

Point très négatif : tout code réalisé est très long, même avec mes nombreuses tentatives d'optimisation.

Prochain mois : 
Initiation aux grammaires et algos associés.

