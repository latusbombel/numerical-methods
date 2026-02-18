reset
set term pngcairo size 800,800 background rgb "white"
set output 'blad2.png'

set view map
set xrange [0:150]
set yrange [0:100]

set xlabel "x"
set ylabel "y"
set title "Błąd (relaksacja lokalna)"
set palette rgbformulae 33,13,10
set colorbox
set size ratio 0.666

plot 'blad.dat' matrix with image