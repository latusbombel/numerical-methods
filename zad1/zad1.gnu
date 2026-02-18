set term png
set out "wynik3_1.png"
set xrange [0:5]
set yrange [-0.1:1]
set ylabel "y(x)"
set xlabel "x"
set title "z.3 - Metoda R. K. 4 y(x)"
plot "wyniki3_1.dat" index 0 u 1:2 w p pt 6 t 'dt = 0.001', "wyniki3_1.dat" index 1 u 1:2 w p pt 1 t 'dt = 0.01', "wyniki3_1.dat" index 2 u 1:2 w p pt 2 t 'dt = 1', "wyniki3_1.dat" index 0 u 1:3 w l lt 1 t 'analityczne'