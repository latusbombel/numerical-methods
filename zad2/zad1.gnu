set term png
set output "wyniki2.png"
set xrange [0:100]
set yrange [0:550]
set ylabel "z(t), u(t)"
set xlabel "t"
set grid
set title "metoda Newtona"
plot 'wyniki2.dat' u 1:2 w p pt 1 t "u(t)", 'wyniki2.dat' u 1:3 w p pt 2 t 'z(t)'