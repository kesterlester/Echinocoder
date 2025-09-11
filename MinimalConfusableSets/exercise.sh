#!/bin/sh

for M in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17
do {

for way in 1 2 3 4
do {

python tom_demo.py $M 6 $way| grep --line-buffered 'AT END.*best_siz'

 echo
} done

echo ----

} done | colrm 100
