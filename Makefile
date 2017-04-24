all: build

build:
	g++ -std=c++11 -Wall src/main.cpp -o main

clean_obj:
	-rm ./main src/*.gch src/*.exe

clean_graph:
	-rm results/*.png

clean_data:
	-rm results/data*.txt

clean_all: clean_obj clean_graph clean_data

complete: clean_obj run plot

pdf:
	-mv docs/report.pdf docs/report_old.pdf
	lyx --export pdf2 docs/report.lyx

plot:
	cd results; gnuplot plot_file.plt;

run: build
	./main
