#!/bin/bash

eigen=$1
cat > plotEnergy.gp << EOF
plot "en1" u 1:2 title "energy" w l, "en" u 1:2 title "FFW" w l, $eigen title "exact ground state energy"
pause 3
reread
EOF
gnuplot plotEnergy.gp
