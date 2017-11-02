#!/bin/bash

eigen=`sort -n eigen.txt | head -n 1`
cat > plotEnergy.gp << EOF
plot "energy.txt" u 1:2 title "energy" w l, "" u 1:3 title "energyMult", $eigen title "exact ground state energy"
pause 3
reread
EOF
#gnuplot plotEnergy.gp
cat > plotCoeff.gp << EOF
#set yrange [0:0.14]
plot "intEigenVec.txt" u 1:(\$2**2), "coeff.txt" u 1:(\$2**2)
pause 1
reread
EOF
gnuplot plotCoeff.gp
