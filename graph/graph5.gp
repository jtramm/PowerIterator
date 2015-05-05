set terminal pdf
set output 'prob5.pdf'
set termoption dash
set style data lines
set key right top Left title 'Legend' box
set title "Normalized Flux of Reactor - Problem 1"
set xlabel "X-coordinate (cm)"
set ylabel "Normalized Flux"
plot 'prob5.dat' using 1:2 title 'Fast' lw 7, 'prob5.dat' using 1:3 title 'Thermal' lw 7 lc 'blue'
