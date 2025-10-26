clear

if [ ! -d "./build" ]; then
	mkdir ./build
fi

if [ -d "./build/src" ]; then
	rm -rf ./build/src
fi

cd build
cmake -B . ..
cmake --build .

cd ./src
./main
