import matplotlib.pyplot as plt
import numpy as np

HT = [0, 0.005, 0.025]   # sec
Q2 = [0.516, 9.185, 197.923]  # sec
h = [0.2, 0.4, 1]  # pas du maillage

plt.plot(h, HT, "o-", label="Table de hachage")
plt.plot(h, Q2, "s-", label="Méthode naïve")

plt.xlabel("Pas du maillage h")
plt.ylabel("Temps d'exécution (s)")
plt.legend()
plt.grid(True)

plt.show()