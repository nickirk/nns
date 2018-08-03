#!/bin/bash

eigen=$1
cat > plotEnergy.gp << EOF
plot "en" u 1:2 title "RBM ADAM 50% random det" w l, "en1" u 1:2 title "RBM SR decaying rand det" w l, "en2" u 1:2 title "RBM ADAM decaying prob to get rand det" w l, "en3" u 1:2  w l title "SR decaying rand det rescaled gradient", $eigen title "exact ground state energy sys 1"
pause 3
reread
EOF
gnuplot plotEnergy.gp
