#!/bin/bash

eigen=`sort eigen.txt | head -n 1`
cat > plot.gp << EOF
plot "energy.txt" u 1:2 title "energy",  $eigen title "exact ground state energy"
pause 1
reread
EOF
gnuplot plot.gp
