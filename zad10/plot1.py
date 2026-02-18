import matplotlib.pyplot as plt
import numpy as np

plik = np.loadtxt("wychylenie3.dat")
t = plik[:,0]
t = np.unique(t)
x = plik[:,1]
x = np.unique(x)
u = plik[:,2]
u = u.reshape(len(t), len(x)).T

X, T = np.meshgrid(x, t)

fig, ax = plt.subplots(figsize=(6, 5))
im = ax.imshow(u, origin="lower", extent=[x.min(), x.max(), t.min(), t.max()], cmap="bwr", aspect="auto")

ax.set_title(r"u(x,t) $\alpha$ = 1, $\beta$ = 1")
ax.set_xlabel("t")
ax.set_ylabel("x")

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("u(x,t)")

plt.savefig("wychylenie3.png")
# print(u[0])