set terminal x11
set size square
set xlabel 'x'
set ylabel 'y'
set xrange [-0.1:1.1]
set yrange [-0.1:1.1]
set view 0,0
splot 'gnuplot.dat'
pause -1