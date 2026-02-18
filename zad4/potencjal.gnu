# set term png
# set output "ped.png"
# # set xrange [0:100]
# # set yrange [0:550]
# # set ylabel "L"
# # set xlabel "t"
# set size square 
# # set grid
# # set title "metoda Newtona"
# plot 'potencjal.dat' u 1:7 w p pt 1
reset
set term pngcairo size 800,800 background rgb "white"
set output 'potencjal2.png'

set view map
set xrange [0:150]
set yrange [0:100]

set xlabel "x"
set ylabel "y"
set title "Potencjał V(x,y) (relaksacja lokalna)"
set palette rgbformulae 33,13,10
set colorbox
set size ratio 0.666

plot 'potencjal.dat' matrix with image