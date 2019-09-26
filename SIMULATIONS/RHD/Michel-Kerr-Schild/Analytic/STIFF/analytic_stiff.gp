#!/bin/bash

gnuplot << EOF

set terminal epslatex color input "" 12 lw 3
set output "${1}.tex"

plot [3:20] sqrt( 1 + 2*(x*x + 2*x + 4)/x**3), \
"${1}.dat" u 1:3

EOF
