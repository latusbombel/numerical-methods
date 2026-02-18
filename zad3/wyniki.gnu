set term png
set output "ped.png"
# set xrange [0:100]
# set yrange [0:550]
set ylabel "L"
set xlabel "t"
set size square 
# set grid
# set title "metoda Newtona"
plot 'wyniki2.dat' u 1:7 w p pt 1