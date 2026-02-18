reset

# --- USTAWIENIA WYJŚCIA ---
set terminal pngcairo size 800,800 background rgb "white"
set output 'potencjal5.png'

# --- OPIS OSI I TYTUŁ ---
set xlabel "x"
set ylabel "y"
set title "Potencjał V(x,y) (siatka k = 1)"

# --- WYGLĄD MAPY ---
set view map
set size ratio -1           # proporcje 1:1
set pm3d map                # rysowanie w trybie mapy 2D
set palette rgbformulae 33,13,10
set colorbox

unset key

# --- WYKRES ---
splot 'potencjal.dat' index 4 using 1:2:3 with pm3d


