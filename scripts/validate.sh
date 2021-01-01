#!/bin/bash
for i in $(seq 0 1 870)
do
   echo "$i"
   ../build/triangulate -s /gisdata/raster/dat/ca_source.dat -v "$i"
done
