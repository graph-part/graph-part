#!/bin/bash

mkdir output

graphpart "$@" >output/signalp6.log.txt 2>&1

#python generate_output.py

#cat output.md