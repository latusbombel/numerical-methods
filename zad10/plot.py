import matplotlib.pyplot as plt
import numpy as np

energia = np.loadtxt("energia3.dat")
# energia1 = np.loadtxt("energia1.dat")
# energia2 = np.loadtxt("energia2.dat")
# energia = np.loadtxt("energia.dat")
t = energia[:,0]
ene1 = energia[:,1]
# ene2 = energia1[:,1]
# ene3 = energia2[:,1]

fig, ax = plt.subplots(figsize = (6,4))
ax.plot(t, ene1, label = r'$\alpha$ = 1, $\beta$ = 1')
# ax.plot(t, ene2, label = r'$\alpha$ = 0, $\beta$ = 0.1')
# ax.plot(t, ene3, label = r'$\alpha$ = 0, $\beta$ = 1')


ax.set_title("E(t)")
ax.set_xlabel("t")
ax.set_ylabel('E(t)')
ax.grid(1)
ax.set_ylim(ene1.min(), ene1.max())

ax.legend(loc="center right", title="Functions", frameon=True, ncol=1)
plt.savefig('energia.png')