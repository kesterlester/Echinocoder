#!/bin/sh
for i in test*.py
do {
  echo Running test $i ==============
  pytest $i
} done
