
# to be executed in the parent directory
cd example
tar -xf p38.tar.xz
cd ..

make

cd example/p38
../../combine.exe -i combine.in -o combine.out

