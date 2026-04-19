# Projet de génération de maillage, MS13 AMS

Tout les codes sources du projet sont dans le fichier **mesh.c**. Les parties 1 et 2 ont été mise en oeuvre dans leur totalité. 
Le compilateur utilisé est MinGw32 (compilateur de windows un peu dépassé), le Makefile à été adapté en conséquence (ne pas hésiter à le remodifier pour faire tourner les codes)

Pour tester les codes, on dispose de 3 fichiers :

- **hsh_table.c** dans lequel il est possible de tester les méthodes implémentées dans la première partie (notamment autour des tables de hachages)
On compile le projet avec la commande suivante : *mingw32-make hsh_table*  et que l'on execute avec *./hsh_table.exe*

-**insertion.c** dans lequel on teste la méthode d'insertion d'un point dans un maillage, que l'on compile avec *mingw32-make insertion* et que l'on exécute avec 
*./insertion.exe*.

- **compression.c** dans lequel on teste les méthodes de compressions d'image, que l'on compile avec *mingw-make compression* et que l'on exécute avec 
*./compression.exe*

Il n'est pas nécéssaire d'ajouter d'arguments dans la commande d'éxécution, les chemins vers les fichiers .mesh sont directement écris dans les scripts. 
Les codes génèrent divers fichier (des .mesh, .sol et .solb), il est recommander entre chaque compilation d'effacer ces fichiers avec la commande 
*mingw32-make clean*.

Si jamais vous rencontrez un soucis lors de la compilation, n'hésitez pas à me contacter par mail. 
