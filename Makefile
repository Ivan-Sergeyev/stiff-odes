all: build_main

build_main:
	g++ -std=c++11 -Wall src/main.cpp -o main

build_test:
	g++ -std=c++11 -Wall src/test.cpp -o test

clean_obj:
	-rm ./main ./test src/*.gch src/*.exe

clean_graph:
	-rm results/*.png

clean_data:
	-rm results/data*.txt

clean_all: clean_obj clean_graph clean_data

complete: clean_obj run plot pdf

pdf:
	-mv docs/report.pdf docs/report_old.pdf
	-lyx --export pdf2 docs/report.lyx

plot:
	cd results; gnuplot plot_file.plt;

run: build_main
	./main

test: build_test
	./test
	cd results_test; gnuplot plot_file.plt;
