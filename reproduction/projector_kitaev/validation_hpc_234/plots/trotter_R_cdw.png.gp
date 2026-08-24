set terminal pngcairo size 900,600
set output '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/validation_hpc_234/plots/trotter_R_cdw.png'
set title 'R_cdw vs dt^2 (correlated jackknife)'
set xlabel 'dt^2'
set ylabel 'R_cdw'
set key outside
f(x)=0.22915104312195334-0.26251274575482419*x
plot '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/validation_hpc_234/plots/trotter_R_cdw.dat' using 1:2:3 with yerrorbars title 'QMC', f(x) title 'correlated GLS fit'
