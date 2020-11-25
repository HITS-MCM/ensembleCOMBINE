
# to be executed in the parent directory
cd example
unzip p38.zip
cd ..

make

cd example/p38
../../combine.exe -i combine.in -o combine.out

