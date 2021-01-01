#!/bin/bash
for i in $(seq 0 1 $2)
do
   echo "$i"
   ../build/triangulate -s $1 -v "$i"
done
