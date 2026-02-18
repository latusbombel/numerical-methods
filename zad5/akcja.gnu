reset
set term png
set output "akcja.png"
# set xrange [0:100]
# set yrange [-100:1000]
set ylabel "S"
set xlabel "it"
# set size square 
# set grid
set title "S(it)"
set logscale x
plot 'akcja.dat' index 0 u 1:2 w p pt 1 title 'k = 16', 'akcja.dat' index 1 u 1:2 w p pt 2 title 'k = 8', 'akcja.dat' index 2 u 1:2 w p pt 3 title 'k = 4', 'akcja.dat' index 3 u 1:2 w p pt 4 title 'k = 2', 'akcja.dat' index 4 u 1:2 w p pt 4 title 'k=1'