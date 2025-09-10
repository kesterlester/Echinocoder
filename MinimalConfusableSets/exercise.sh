#!/bin/sh

for M in 1 2 3 4 5 6 7 8 9 10; do { python tom_demo.py $M 2 | grep --line-buffered 'AT END.*smallest' ; echo ; } done | colrm 100
