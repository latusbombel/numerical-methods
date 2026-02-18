import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# --- 1. KONFIGURACJA I WCZYTANIE SIATKI ---
# Wczytujemy jeden plik tylko po to, żeby ustalić wymiary (nx, ny) i osie (X, Y)
# Zakładam, że plik0.dat istnieje. Jeśli nie, użyj predkosc.dat do ustalenia siatki.
data_init = np.loadtxt('predkosc.dat') 

# Zakładamy format: x y wartość
x_flat = data_init[:, 0]
y_flat = data_init[:, 1]

unique_x = np.unique(x_flat)
unique_y = np.unique(y_flat)

nx = len(unique_x)
ny = len(unique_y)

# Tworzenie siatki
X, Y = np.meshgrid(unique_x, unique_y)

# --- 2. PRZYGOTOWANIE WYKRESU ---
fig, ax = plt.subplots(figsize=(10, 4))

# Inicjalizacja pustą macierzą zer lub pierwszym plikiem
# vmin i vmax są KLUCZOWE, żeby kolory nie "mrugały" przy zmianie wartości
# Dostosuj vmax do maksymalnej spodziewanej wartości (np. 1.0 dla gęstości)
im = ax.pcolormesh(X, Y, np.zeros((ny, nx)), cmap="plasma", shading='auto', vmin=0, vmax=1.0)

ax.set_title("Symulacja Adwekcji-Dyfuzji")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect('auto') # Rozciąga wykres na całe okno

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("Gęstość u")

# --- 3. FUNKCJA AKTUALIZUJĄCA ---
def animate(i):
    # 1. Konstrukcja nazwy pliku (zgodna z Twoim kodem C++)
    filename = f"gify/gestosc{i}.dat" 
    
    try:
        # 2. Wczytanie danych dla bieżącej klatki
        data = np.loadtxt(filename)
        vals_flat = data[:, 2] # Zakładam, że gęstość jest w 3 kolumnie
        
        # 3. Reshape i transpozycja (tak samo jak robiłeś wcześniej)
        vals_2d = vals_flat.reshape(nx, ny).T
        
        # 4. Aktualizacja wykresu
        # UWAGA: pcolormesh.set_array oczekuje spłaszczonej tablicy (ravel) 
        # i to bez wymiarów X, Y (same wartości wewnątrz)
        im.set_array(vals_2d.ravel())
        
        ax.set_title(f"Krok czasowy: {i}")
        
    except OSError:
        print(f"Nie znaleziono pliku: {filename}")
        # Można tu przerwać animację lub pominąć klatkę
    
    return [im]

# --- 4. TWORZENIE I ZAPIS ANIMACJI ---
# frames=50 -> ile masz plików (np. od plik0.dat do plik49.dat)
n_frames = 50 

ani = animation.FuncAnimation(fig, animate, frames=n_frames, interval=100, blit=False)

print("Generowanie GIFa...")
writer = animation.PillowWriter(fps=10, metadata=dict(artist='Me'), bitrate=1800)
ani.save('symulacja.gif', writer=writer)
print("Gotowe!")