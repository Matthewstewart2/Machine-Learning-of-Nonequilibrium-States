#/bin/bash

#2'000 steps printing on screen
# time and magnetisation
#Phase transition at T=0

temp=1;

c++ Ising1d.c++ -o Ising1d -O3
./Ising1d ${temp}

gnuplot <<- EOF
L=100
load "plot.plt"
EOF


rm out*
