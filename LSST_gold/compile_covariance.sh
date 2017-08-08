#!/bin/bash 

gcc -O3 -std=c99 -Wno-missing-braces -Wno-missing-field-initializers -I$HOME/local/include -L$HOME/local/lib  -o ./compute_covariances_fourier compute_covariances_fourier.c -lfftw3 -lgsl -lgslcblas -lm -L../../cosmolike_core/class -lclass
