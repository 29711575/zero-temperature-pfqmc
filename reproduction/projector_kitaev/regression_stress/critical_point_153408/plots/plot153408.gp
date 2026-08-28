set terminal pngcairo size 1100,750 enhanced font 'Arial,16'
set xlabel 'V'; set ylabel 'R_cdw'; set grid; set key left top
set title '153408 known-physics: R_cdw(V), OBC hs=1, theta=L'
set output 'R_cdw_crossings.png'
plot 'R_cdw_plot.dat' index 0 using 2:3:4 with yerrorlines lw 2 pt 7 title 'L=10', \
     '' index 1 using 2:3:4 with yerrorlines lw 2 pt 7 title 'L=18', \
     '' index 2 using 2:3:4 with yerrorlines lw 2 pt 7 title 'L=26', \
     '' index 3 using 2:3:4 with yerrorlines lw 2 pt 7 title 'L=34'
set xlabel '(V-V_c)L'; set ylabel 'R_cdw'; set key left top
set title '153408 data collapse (fixed nu=1, central fit)'
set output 'R_cdw_collapse.png'
plot 'collapse_plot.dat' using 1:2:3 with yerrorbars pt 7 title 'theta=L'
