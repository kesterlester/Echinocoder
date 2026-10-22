#!/bin/sh

for i in *list_of_reals*.py ; do { echo ====== ; echo $i ; echo ====== ; python $i ; } done


