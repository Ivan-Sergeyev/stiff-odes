all: build

build:
	g++ -std=c++11 src/main.cpp -o main

clean:
	-rm ./main src/*.gch results/*.png

complete: clean run plot

plot:
	cd results; gnuplot plot_file.plt;

run: build
	./main
