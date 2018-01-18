#!/bin/bash

eigen=-5.95
cat > plotEnergy.gp << EOF
plot "energy_gd.txt" u 1:2 title "energy 20" w l, "energy_test.txt" u 1:2 title "test", "energy_40hidden.txt" u 1:2 w l title "energy hidden 40", "energy_RMSprop.txt" u 1:2 w l, $eigen title "exact ground state energy"
pause 3
reread
EOF
gnuplot plotEnergy.gp
