#!/bin/sh

for M in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
do {

for way in 4 #1 2 3 4 5
do {

python tom_demo.py $M 2 -o | grep --line-buffered 'AT END.*\(\(best_siz\)\|\(config\)\)'

 echo
} done

echo ----

} done | colrm 100
