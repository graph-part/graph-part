#!/bin/bash

mkdir output

graphpart "$@" >output/graphpart.log.txt #2>&1
python make_output.py
#python generate_output.py

cat output.md