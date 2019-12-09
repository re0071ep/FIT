set terminal windows 0 color solid enhanced font "Meiryo, 9" size 600,500


set linestyle 1 lt 1 lw 1 lc rgb "navy" pt 7 ps 0.6
set linestyle 2 lt 1 lw 2 lc rgb "dark-pink" pt 7 ps 0.4
set linestyle 3 lt 1 lw 1 lc rgb "green" pt 7 ps 0.6
set linestyle 4 lt 1 lw 1 lc rgb "red" pt 7 ps 0.6
set linestyle 5 lt 1 lw 1 lc rgb "blue" pt 7 ps 0.6
set linestyle 6 lt 1 lw 1 lc rgb "orange" pt 7 ps 0.6
set linestyle 7 lt 1 lw 1 lc rgb "dark-green" pt 7 ps 0.6
set linestyle 8 lt 1 lw 1 lc rgb "light-blue" pt 7 ps 0.6
set linestyle 9 lt 1 lw 1 lc rgb "dark-red" pt 7 ps 0.6
set linestyle 10 lt 1 lw 1 lc rgb "dark-olivegreen" pt 7 ps 0.6

set linestyle 11 lt 1 lw 1 lc rgb "navy" pt 2 ps 0.6
set linestyle 12 lt 1 lw 1 lc rgb "dark-pink" pt 2 ps 0.6
set linestyle 13 lt 1 lw 1 lc rgb "green" pt 2 ps 0.6
set linestyle 14 lt 1 lw 1 lc rgb "red" pt 2 ps 0.6
set linestyle 15 lt 1 lw 1 lc rgb "blue" pt 2 ps 0.6
set linestyle 16 lt 1 lw 1 lc rgb "orange" pt 2 ps 0.6
set linestyle 17 lt 1 lw 1 lc rgb "dark-green" pt 2 ps 0.6
set linestyle 18 lt 1 lw 1 lc rgb "light-blue" pt 2 ps 0.6
set linestyle 19 lt 1 lw 1 lc rgb "dark-red" pt 2 ps 0.6
set linestyle 20 lt 1 lw 1 lc rgb "dark-olivegreen" pt 2 ps 0.6

set linestyle 21 lt 1 lw 1 lc rgb "navy" pt 10 ps 0.6
set linestyle 22 lt 1 lw 1 lc rgb "dark-pink" pt 10 ps 0.6
set linestyle 23 lt 1 lw 1 lc rgb "green" pt 10 ps 0.6
set linestyle 24 lt 1 lw 1 lc rgb "red" pt 10 ps 0.6
set linestyle 25 lt 1 lw 1 lc rgb "blue" pt 10 ps 0.6
set linestyle 26 lt 1 lw 1 lc rgb "orange" pt 10 ps 0.6
set linestyle 27 lt 1 lw 1 lc rgb "dark-green" pt 10 ps 0.6
set linestyle 28 lt 1 lw 1 lc rgb "light-blue" pt 10 ps 0.6
set linestyle 29 lt 1 lw 1 lc rgb "dark-red" pt 10 ps 0.6
set linestyle 30 lt 1 lw 1 lc rgb "dark-olivegreen" pt 10 ps 0.6

#set lmargin 0
#set rmargin 0
#set tmargin 1
#set bmargin 3

set grid

##### plot1 ####
set xrange [ *: * ] noreverse nowriteback
set yrange [ * : * ] noreverse nowriteback
set logscale x
set xtics format ""
set ytics format ""
set xtics add ("10^{-2}" 1e-2,"0.1" 1e-1,"1" 1,"10" 1e1,"10^{2}" 1e2,"10^{3}" 1e3,"10^{4}" 1e4,"10^{5}" 1e5,"10^{6}" 1e6,"10^{7}" 1e7)
set ytics add ("0.1" 1e-1,"1" 1,"10" 1e1,"10^{2}" 1e2,"10^{3}" 1e3,"10^{4}" 1e4,"10^{5}" 1e5,"10^{6}" 1e6,"10^{7}" 1e7,"10^{8}" 1e8,"10^{9}" 1e9,"10^{10}" 1e10,"10^{11}" 1e11)
set xlabel "f [Hz]" 
set ylabel "|Z| [{/Symbol W}]" 
set key right top
set size ratio 0.65
set datafile separator ","

plot \
"gnuplot1.csv" using 1:(-$3) title "Experiment 1" with point ls 1 \
,"gnuplot2.csv" using 1:(-$3) title "Experiment 2" with point ls 2 \
,"gnuplot3.csv" using 1:(-$3) title "Experiment 3" with point ls 3 \
,"gnuplot4.csv" using 1:(-$3) title "Experiment 4" with point ls 4 \
,"gnuplot5.csv" using 1:(-$3) title "Experiment 5" with point ls 5 \
,"gnuplot6.csv" using 1:(-$3) title "Experiment 6" with point ls 6 \
,"gnuplot7.csv" using 1:(-$3) title "Experiment 7" with point ls 7 \
,"gnuplot8.csv" using 1:(-$3) title "Experiment 8" with point ls 8 \
,"gnuplot9.csv" using 1:(-$3) title "Experiment 9" with point ls 9 \
,"gnuplot10.csv" using 1:(-$3) title "Experiment 10" with point ls 10 \
,"gnuplot1.csv" using 1:(-$5) title "Fitting 1" with point ls 11 \
,"gnuplot2.csv" using 1:(-$5) title "Fitting 2" with point ls 12 \
,"gnuplot3.csv" using 1:(-$5) title "Fitting 3" with point ls 13 \
,"gnuplot4.csv" using 1:(-$5) title "Fitting 4" with point ls 14 \
,"gnuplot5.csv" using 1:(-$5) title "Fitting 5" with point ls 15 \
,"gnuplot6.csv" using 1:(-$5) title "Fitting 6" with point ls 16 \
,"gnuplot7.csv" using 1:(-$5) title "Fitting 7" with point ls 17 \
,"gnuplot8.csv" using 1:(-$5) title "Fitting 8" with point ls 18 \
,"gnuplot9.csv" using 1:(-$5) title "Fitting 9" with point ls 19 \
,"gnuplot10.csv" using 1:(-$5) title "Fitting 10" with point ls 20