#!/bin/bash

eigen=`sort -n eigen.txt | head -n 1`
cat > plotEnergy.gp << EOF
plot "energy.txt" u 1:2 title "energy" w l,  -3.66871 title "exact ground state energy"
pause 3
reread
EOF
gnuplot plotEnergy.gp
