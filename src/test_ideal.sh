#!/bin/bash
DATA_PATH="/home/qmh/data"
BUILD_PATH="/home/qmh/IDEAL/build"
DATA1="has_child.idl"
DATA2="child.idl"
DATA3="sampled.points.dat"
> output.txt
make clean
make contain -j16
$BUILD_PATH/contain -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA3 -r >> output.txt
$BUILD_PATH/contain -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA3 --vector >> output.txt
$BUILD_PATH/contain -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA3 -q >> output.txt
make clean
make contain_polygon -j16
$BUILD_PATH/contain_polygon -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA2 -r >> output.txt
$BUILD_PATH/contain_polygon -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA2 --vector >> output.txt
$BUILD_PATH/contain_polygon -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA2 -q >> output.txt
make clean
make within -j16
$BUILD_PATH/within -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA3 -r >> output.txt
$BUILD_PATH/within -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA3 --vector >> output.txt
$BUILD_PATH/within -s $DATA_PATH/$DATA1 -t $DATA_PATH/$DATA3 -q >> output.txt
make clean
make within_polygon -j16
$BUILD_PATH/within_polygon -s $DATA_PATH/$DATA1 -r >> output.txt
$BUILD_PATH/within_polygon -s $DATA_PATH/$DATA1 --vector >> output.txt
$BUILD_PATH/within_polygon -s $DATA_PATH/$DATA1 -q >> output.txt



