make clean
make within USE_GPU=1 -j16
# ../build/within -s /home/qmh/data/complex.idl -t /home/qmh/data/points.dat -r
../build/within -s /home/qmh/data/has_child.idl -t /home/qmh/data/child_points.dat -r
