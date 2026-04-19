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

## Dummy comp

RCT_dummy = [0.011898, 0.142972, 0.333636, 0.508873]
time_dummy = [0.010000, 0.175000, 0.582000, 1.259000]

RCT_random = [0.100020, 0.295120, 0.490623, 0.699536]
time_random = [0.085000, 0.514000, 1.216000, 1.875000]

RCT_psnr = [0.003428, 0.012704, 0.031760, 0.048599 ]
time_psnr = [0.496000, 3.513000, 13.398000, 23.775000 ]

plt.plot(RCT_dummy, time_dummy, label = 'Dummy Compresser')
plt.plot(RCT_random, time_random, label = 'Random Compresser')
plt.plot(RCT_psnr, time_psnr, label = 'PNSR Compresser')

plt.xlabel('Ratio de compression')
plt.ylabel("Temps d'exécution")
plt.grid()
plt.legend()
plt.show()


PNSR_nat = [20.228425, 25.388909, 30.039976, 35.004725, 38.009816]
time_nat = [0.724000, 7.716000,23.056000, 66.692000, 133.438000]

PNSR_opt = [20.228425, 25.388909, 30.422689]
time_opt = [0.612000, 6.130000, 21.121000]

plt.plot(PNSR_nat, time_nat, label = 'Algorithme non optimisé')
plt.plot(PNSR_opt, time_opt, label = 'Algorithme optimisé')

plt.title('Algo optimisé contre non optimisé')
plt.xlabel("PSNR")
plt.ylabel("Teps d'exécution")
plt.legend()
plt.grid()

plt.show()