reset
set term gif size 800,300 animate delay 10
set output "anim.gif"
n=49    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]
# set dgrid3d 400,90

do for [i=0:n] {
  file = sprintf("gify/gestosc%i.dat",i)
  splot file u 1:2:3 w pm3d  title sprintf("t=%i",i)
} 