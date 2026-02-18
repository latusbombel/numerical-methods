reset
set term png
set output "akcja2.png"
# set xrange [0:100]
# set yrange [-100:1000]
set ylabel "S"
set xlabel "it"
# set size square 
# set grid
set title "S(it) (relaksacja lokalna)"
set logscale x
plot 'akcja.dat' index 0 u 1:2 w p pt 1 title '{/symbol w} = 1.0', 'akcja.dat' index 1 u 1:2 w p pt 2 title '{/symbol w} = 1.4', 'akcja.dat' index 2 u 1:2 w p pt 3 title '{/symbol w} = 1.8', 'akcja.dat' index 3 u 1:2 w p pt 4 title '{/symbol w} = 1.9'