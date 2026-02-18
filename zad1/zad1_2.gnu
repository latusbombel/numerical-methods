set term png
set out "wynik3_2.png"
set yrange [-0.05:0.05]
set xrange [0:5]
set ylabel "błąd(x)"
set xlabel "x"
set title "z.3 - Metoda R. K. 4 błąd (x)"
plot "wyniki3_2.dat" index 0 u 1:2 w p pt 1 t "błąd dt = 0.01", "wyniki3_2.dat" index 1 u 1:2 w p pt 2 t "błąd dt = 0.1", "wyniki3_2.dat" index 2 u 1:2 w p pt 3 t "błąd dt = 1"