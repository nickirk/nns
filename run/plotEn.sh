#!/bin/bash

eigen=$1
cat > plotEnergy.gp << EOF
plot "en" u 1:2 title "FFW FullSam" w l, "en1" u 1:2 title "Ab" w l, "en2" u 1:2 title "RBM fullSam" w l, $eigen title "exact ground state energy sys 1"
pause 3
reread
EOF
gnuplot plotEnergy.gp
