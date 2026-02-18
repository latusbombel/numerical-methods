import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('calki.dat')
data2 = np.loadtxt('calki1.dat')

t = data[:, 0]
xSr = data[:, 1]
c = data[:, 2]

xSr1 = data2[:, 1]
c1 = data2[:, 2]

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(t, xSr, label = "$x_{sr}$, D = 0")
ax.plot(t, c, label = "c, D = 0")
ax.plot(t, c1,label = "c, D = 0.1")
ax.plot(t, xSr1, label = "$x_{sr}$, D = 0.1")

ax.set_title("Wykres x średniego i całki gęstości")
ax.set_xlabel("t")
ax.set_ylabel('$x_{sr}$, c(t)')
ax.grid(1)

ax.legend(loc="upper right", title="Functions", frameon=True, ncol=1)
plt.savefig('wynik2.png')