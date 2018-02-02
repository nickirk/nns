#!/bin/bash

eigen=-2.1027
cat > plotEnergy.gp << EOF
plot "energy" u 1:2 title "energy" w l,  $eigen title "exact ground state energy"
pause 3
reread
EOF
gnuplot plotEnergy.gp
