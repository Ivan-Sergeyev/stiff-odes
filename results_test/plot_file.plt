reset

set terminal pngcairo dashed enhanced font "Sans,15" size 1000,600
set encoding utf8
set grid
set key top left

# todo: add titles?

set style line 1 lt 1 lc 1 lw 4
set style line 2 lt 1 lc 2 lw 3
set style line 3 lt 1 lc 3 lw 2
set style line 4 lt 3 lc 4 lw 1

datafile1 = "data_1.txt"
datafile2 = "data_2.txt"
datafile3 = "data_3.txt"

ttl1 = "rk1 expl"
ttl2 = "rk4 expl"
ttl3 = "rk4 impl"

do for [col=2:5] {
    set output "graph_test_x".(col - 1)."_t.png"
    set ylabel "x_".(col - 1) offset 2,0,0
    set xlabel "t" offset 0,0.5,0

    p = 2**(col - 4)

    plot datafile1 using 1:col with lines ls 1 title ttl1,\
         datafile2 using 1:col with lines ls 2 title ttl2,\
         datafile3 using 1:col with lines ls 3 title ttl3,\
         p*exp(x)              with lines ls 4 title sprintf("%0.1f*exp(x)", p)
}
