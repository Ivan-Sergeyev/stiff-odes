reset

set terminal pngcairo dashed enhanced font "Sans,11" size 1000,600
set encoding utf8
set grid
set key top left

# todo: add titles

do for [m=1:2] {  # m=1:3
    do for [i=1:3] {
    	datafile = "data_".i."_".m.".txt"

    	set output "graph_".i."_".m."_x.png"
    	set xlabel "t" offset 0,0.5,0
    	set ylabel "x" offset 2,0,0
    	plot datafile using 1:2 with points lw 2 lc 3 lt 1 notitle

    	set output "graph_".i."_".m."_y.png"
    	set xlabel "t" offset 0,0.5,0
    	set ylabel "y" offset 2,0,0
    	plot datafile using 1:3 with points lw 2 lc 3 lt 1 notitle

        set output "graph_".i."_".m."_xy.png"
        set xlabel "x" offset 0,0.5,0
        set ylabel "y" offset 2,0,0
        plot datafile using 2:3 with points lw 2 lc 3 lt 1 notitle
    }

    datafile = "data_1_".m.".txt"

    set output "graph_1_".m."_a.png"
    set xlabel "t" offset 0,0.5,0
    set ylabel "α" offset 2,0,0
    plot datafile using 1:4 with points lw 2 lc 3 lt 1 notitle

    do for [i=2:3] {
        datafile = "data_".i."_".m.".txt"

        set output "graph_".i."_".m."_a1.png"
        set xlabel "t" offset 0,0.5,0
        set ylabel "α_1" offset 2,0,0
        plot datafile using 1:4 with points lw 2 lc 3 lt 1 notitle

        set output "graph_".i."_".m."_a2.png"
        set xlabel "t" offset 0,0.5,0
        set ylabel "α_2" offset 2,0,0
        plot datafile using 1:5 with points lw 2 lc 3 lt 1 notitle
    }
}
