reset

set encoding utf8
set grid
set key top left

# todo: add titles?

col_t = 1
col_x = 2
col_y = 3
col_a = 4
col_a1 = 4
col_a2 = 5

set style line 1 lt 1 lc 3 lw 2

do for [i=1:3] {
    do for [m=1:3] {
    	datafile = "data_".i."_".m.".txt"

        if (i == 1) {
            set terminal pngcairo dashed enhanced font "Sans,15" size 2800,750
        } else {
            set terminal pngcairo dashed enhanced font "Sans,15" size 2800,600
        }

    	set output "graph_".i."_".m."_x_t.png"
        set ylabel "x" offset 2,0,0
    	set xlabel "t" offset 0,0.5,0

        plot datafile every 2 using col_t:col_x with lines ls 1 notitle

    	set output "graph_".i."_".m."_y_t.png"
        set ylabel "y" offset 2,0,0
    	set xlabel "t" offset 0,0.5,0
        plot datafile every 2 using col_t:col_y with lines ls 1 notitle

        if (i == 1) {
            set output "graph_".i."_".m."_a_t.png"
            set ylabel "α" offset 2,0,0
            set xlabel "t" offset 0,0.5,0
            plot datafile every 2 using col_t:col_a with lines ls 1 notitle
        } else {
            set output "graph_".i."_".m."_a1_t.png"
            set ylabel "α_1" offset 2,0,0
            set xlabel "t" offset 0,0.5,0
            plot datafile every 2 using col_t:col_a1 with lines ls 1 notitle

            set output "graph_".i."_".m."_a2_t.png"
            set ylabel "α_2" offset 2,0,0
            set xlabel "t" offset 0,0.5,0
            plot datafile every 2 using col_t:col_a2 with lines ls 1 notitle
        }

        if (i == 1) {
            set terminal pngcairo dashed enhanced font "Sans,15" size 750,750
        } else {
            set terminal pngcairo dashed enhanced font "Sans,15" size 600,600
        }

        set output "graph_".i."_".m."_y_x.png"
        set ylabel "y" offset 2,0,0
        set xlabel "x" offset 0,0.5,0
        plot datafile every 2 using col_x:col_y with lines ls 1 notitle

        if (i == 1) {
            set output "graph_".i."_".m."_x_a.png"
            set ylabel "x" offset 2,0,0
            set xlabel "α" offset 0,0.5,0
            plot datafile every 2 using col_a:col_x with lines ls 1 notitle

            set output "graph_".i."_".m."_y_a.png"
            set ylabel "y" offset 2,0,0
            set xlabel "α" offset 0,0.5,0
            plot datafile every 2 using col_a:col_y with lines ls 1 notitle
        } else {
            set output "graph_".i."_".m."_x_a1.png"
            set ylabel "x" offset 2,0,0
            set xlabel "α_1" offset 0,0.5,0
            plot datafile every 2 using col_a1:col_x with lines ls 1 notitle

            set output "graph_".i."_".m."_y_a2.png"
            set ylabel "y" offset 2,0,0
            set xlabel "α_2" offset 0,0.5,0
            plot datafile every 2 using col_a2:col_y with lines ls 1 notitle

            set output "graph_".i."_".m."_a2_a1.png"
            set ylabel "α_2" offset 2,0,0
            set xlabel "α_1" offset 0,0.5,0
            plot datafile every 2 using col_a1:col_a2 with lines ls 1 notitle
        }
    }
}
