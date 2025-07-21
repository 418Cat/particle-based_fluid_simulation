clear
rm -rf build/src
cd build
cmake -B . ..
cmake --build .

cd ./src
./main
