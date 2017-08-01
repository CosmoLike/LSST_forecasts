#!/bin/bash 

rm ./like_fourier

gcc -O3 -Wno-missing-braces -Wno-missing-field-initializers -I$HOME/local/include -L$HOME/local/lib  -o ./like_fourier like_fourier.c -lfftw3 -lgsl -lgslcblas -lm -L../../cosmolike_core/class -lclass  -std=gnu99

./like_fourier


