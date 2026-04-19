import matplotlib.pyplot as plt


x = [20, 40, 60, 80, 100, 120, 140]
y = [0.043000, 0.120000, 0.128000, 0.166000, 0.211000, 0.237000, 0.316000]

plt.plot(x,y, label = "Performance CPU")
plt.xlabel("nombre de points insérés")
plt.ylabel("Temps CPU")
plt.legend()
plt.title("Performance CPU de l'insertion de point")
plt.grid(True)
plt.show()

