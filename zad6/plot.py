import matplotlib.pyplot as plt
import numpy as np

# 1. WCZYTANIE DANYCH
# Zakładam, że plik nazywa się 'wyniki.txt' i ma format: numer x y wartość
# Jeśli nie masz pliku, sekcja poniżej wygeneruje przykładowe dane, żeby kod działał.
try:
    data = np.loadtxt('potencjal.dat')
except OSError:
    print("Nie znaleziono pliku 'wyniki.txt'")

# 2. PRZYGOTOWANIE DANYCH DO WYKRESU
# Rozpakowanie kolumn
x_raw = data[:, 1]
y_raw = data[:, 2]
z_raw = data[:, 3]

# Wyznaczenie wymiarów siatki (nx, ny)
unique_x = np.unique(x_raw)
unique_y = np.unique(y_raw)
nx = len(unique_x)
ny = len(unique_y)

# Przekształcenie płaskich wektorów w macierze 2D
# Uwaga: reshape zależy od kolejności zapisu w pętli (najpierw po x czy po y).
# Jeśli obrazek wyjdzie obrócony, zamień 'order' na 'F' lub 'C'.
X = x_raw.reshape((ny, nx))
Y = y_raw.reshape((ny, nx))
Z = z_raw.reshape((ny, nx))

# 3. RYSOWANIE WYKRESU
fig, ax = plt.subplots(figsize=(6, 5))

# Używamy pcolormesh (najlepszy do map kolorów na siatkach)
# cmap='bwr' to Blue-White-Red (taka jak na obrazku)
# vmin=-10, vmax=10 ustawia sztywny zakres kolorów (niebieski=-10, czerwony=10)
mesh = ax.pcolormesh(X, Y, Z, cmap='bwr', shading='auto', vmin=-1, vmax=1)

# 4. DODATKI (Pasek kolorów, Tytuły, Osie)
cbar = fig.colorbar(mesh, ax=ax)
cbar.set_label('V', fontsize=12)
# Ustawienie ticków na colorbarze, żeby było widać 0 pośrodku
cbar.set_ticks([-1, 0, 1])

ax.set_title(f'nx=ny=100, $\epsilon_1=1$, $\epsilon_2=10$', fontsize=12)
ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)

# Ustawienie "ciasnych" granic osi, żeby wykres dotykał ramki
ax.set_xlim(np.min(x_raw), np.max(x_raw))
ax.set_ylim(np.min(y_raw), np.max(y_raw))

# Ustawienie proporcji osi na równe (żeby kwadrat był kwadratem)
ax.set_aspect('equal')

# Opcjonalnie: Zwiększenie czcionki na osiach
ax.tick_params(axis='both', which='major', labelsize=12)

plt.tight_layout()
plt.show()