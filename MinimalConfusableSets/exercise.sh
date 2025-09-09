#!/bin/sh

for M in 1 2 3 4 5 6 7; do { python tom_demo.py $M 2 | grep 'AT END' ; } done | colrm 100
