#!/bin/sh

for i in *array_of_reals*.py ; do { echo ====== ; echo $i ; echo ====== ; python $i ; } done


