#!/bin/bash

eigen=$1
cat > plotEnergy.gp << EOF
plot "en" u 1:2 title "energy" w l,  $eigen title "exact ground state energy"
pause 3
reread
EOF
gnuplot plotEnergy.gp
