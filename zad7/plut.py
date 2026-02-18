import numpy as np
import matplotlib.pyplot as plt

def plot_fluid_data(filename):
    try:
        # 1. Wczytywanie danych z pliku tekstowego
        # Loadtxt automatycznie obsługuje spacje i nowe linie
        data = np.loadtxt(filename)
    except OSError:
        print(f"Błąd: Nie znaleziono pliku '{filename}'. Upewnij się, że jest w tym samym folderze.")
        return

    # 2. Pobranie wymiarów siatki na podstawie danych
    ny, nx = data.shape
    
    # Parametr delta z Twojego kodu C++ (potrzebny do skalowania osi)
    delta = 0.01
    
    # Tworzenie siatki współrzędnych (żeby osie miały fizyczne wartości w metrach, a nie indeksy)
    x = np.linspace(0, (nx - 1) * delta, nx)
    y = np.linspace(0, (ny - 1) * delta, ny)
    X, Y = np.meshgrid(x, y)

    # 3. Rysowanie wykresu
    plt.figure(figsize=(10, 6))
    
    # Wykres konturowy wypełniony (contourf)
    # levels=50 określa gładkość przejść kolorów
    # cmap='jet' lub 'viridis' to schematy kolorów
    contour = plt.contourf(X, Y, data, levels=60, cmap='jet')
    
    # Dodanie paska legendy kolorów
    cbar = plt.colorbar(contour)
    cbar.set_label('Wartość')

    # Ustawienia osi
    plt.title(f"Wizualizacja danych z pliku {filename}")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    
    # Ważne: origin='lower' w imshow/contourf jest domyślny dla X,Y meshgrid,
    # ale upewniamy się, że proporcje są 1:1 (żeby kanał nie był rozciągnięty)
    plt.axis('scaled')

    # Opcjonalnie: Rysowanie szarego prostokąta w miejscu przeszkody
    # Zgodnie z C++: i1=50, j1=55 => x=0.5, y=0.55
    obstacle = plt.Rectangle((0, 0), 50*delta, 55*delta, facecolor='gray', edgecolor='black')
    plt.gca().add_patch(obstacle)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Tutaj wpisz nazwę swojego pliku z danymi
    plot_fluid_data("psi.dat")