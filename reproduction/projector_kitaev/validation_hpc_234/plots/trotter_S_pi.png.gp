set terminal pngcairo size 900,600
set output '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/validation_hpc_234/plots/trotter_S_pi.png'
set title 'S_pi vs dt^2 (correlated jackknife)'
set xlabel 'dt^2'
set ylabel 'S_pi'
set key outside
f(x)=0.045795058937424887-0.020978338888210768*x
plot '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/validation_hpc_234/plots/trotter_S_pi.dat' using 1:2:3 with yerrorbars title 'QMC', f(x) title 'correlated GLS fit'
