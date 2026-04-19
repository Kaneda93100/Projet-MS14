import matplotlib.pyplot as plt

### psnr optimised
iter = [50, 100, 150, 200, 250, 300, 350, 400]

perf_loop_opt = [0.18, 0.24, 0.324, 0.329, 0.32, 0.457, 0.464, 0.507]
perf_pnsr_opt = [0.08, 0.124, 0.173, 0.2150, 0.191, 0.231, 0.242, 0.245]

perf_loop = [0.249, 0.266, 0.352, 0.37, 0.381, 0.424, 0.565, 0.562]
perf_psnr = [0.253, 0.24, 0.347, 0.352, 0.338, 0.341, 0.436, 0.644]

fig1 = plt.figure(figsize=(10,10))
plt.plot(iter, perf_loop, label = 'No opt')
plt.plot(iter, perf_loop_opt, label = 'Opt')
plt.title("Boucle de calcul de l'EQM max")
plt.xlabel("itération")
plt.ylabel("Temps d'exécution")
plt.grid()
plt.legend()

fig2 = plt.figure(figsize=(10,10))
plt.plot(iter, perf_psnr, label = 'No opt')
plt.plot(iter, perf_pnsr_opt, label = 'Opt')
plt.title("Calcul du PNSR")
plt.xlabel("itération")
plt.ylabel("Temps d'exécution")
plt.legend()
plt.grid()

plt.show()


