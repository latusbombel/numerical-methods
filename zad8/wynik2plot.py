import numpy as np
import matplotlib.pyplot as plt

# 1. Wczytanie danych
data = np.loadtxt('predkosc.dat')

# Wyciągnięcie kolumn (to są tablice 1D)
x_flat = data[:, 0]
y_flat = data[:, 1]
vx_flat = data[:, 3]

# 2. Wyznaczenie wymiarów siatki i unikalnych osi
unique_x = np.unique(x_flat)
unique_y = np.unique(y_flat)

nx = len(unique_x)
ny = len(unique_y)

print(f"Wykryto siatkę: {nx} x {ny}") # Powinno wypisać np. 400 x 90

# 3. Zmiana kształtu danych (Reshape)
# UWAGA: Kolejność reshape zależy od tego, jak zapisywałeś pętle w C++.
# Jeśli pętla była: for i (x)... for j (y)..., to dane są ułożone blokami po Y.
# reshape(nx, ny) tworzy macierz, gdzie wiersze to X, a kolumny to Y.
# Musimy ją transponować (.T), żeby Y było na osi pionowej (wiersze), a X na poziomej (kolumny).

vx_2d = vx_flat.reshape(nx, ny).T 

# Tworzenie poprawnego meshgridu (z unikalnych wartości!)
X, Y = np.meshgrid(unique_x, unique_y)

# 4. Rysowanie
fig, ax = plt.subplots(figsize=(10, 4)) # Rozciągnięty wykres, bo kanał jest długi

# shading='auto' lub 'nearest' jest bezpieczniejsze przy pcolormesh
im = ax.pcolormesh(X, Y, vx_2d, cmap="plasma", shading='auto')

ax.set_title("Mapa prędkości $v_y$")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect('auto') # Ważne, żeby zachować proporcje geometryczne

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("$v_y$")

plt.savefig('vy.png')
# plt.show()