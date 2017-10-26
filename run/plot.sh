#!/bin/bash

eigen=`sort -n eigen.txt | head -n 1`
cat > plot.gp << EOF
plot "energy.txt" u 1:2 title "energy" w l, "" u 1:3 title "sampleEnergy" w l, $eigen title "exact ground state energy"
pause 3
reread
EOF
gnuplot plot.gp
