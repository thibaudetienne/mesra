#!/bin/bash

defaultx=`echo 0.50 | bc -l`
x=${1:-$defaultx}

echo $x
