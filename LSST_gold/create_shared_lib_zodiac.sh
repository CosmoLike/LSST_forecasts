#!/bin/bash 

rm like_fourier.so

gcc -Wno-missing-braces -Wno-missing-field-initializers  -I$HOME/local/include -L$HOME/local/lib  -shared -o like_fourier.so -fPIC like_fourier.c -lfftw3 -lgsl -lgslcblas -lm  -L../../cosmolike_core/class -lclass  -std=gnu99
