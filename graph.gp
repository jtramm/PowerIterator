set style data lines
set key right top Left title 'Legend' box
set title "Flux in Reactor"
set xlabel "X-coordinate (cm)"
set ylabel "Normalized Flux"
plot 'data.dat' using 1:2 title 'Fast', 'data.dat' using 1:3 title 'Thermal'
