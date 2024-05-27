make clean
make contain_polygon USE_GPU=1 -j16
../build/contain_polygon -s /home/qmh/data/has_child.idl -t /home/qmh/data/child.idl -r -n 1
