#!/bin/bash

eigen=-4.6984
cat > plotEnergy.gp << EOF
plot "energy.txt" u 1:2 title "energy" w l,  $eigen title "exact ground state energy"
pause 3
reread
EOF
gnuplot plotEnergy.gp
