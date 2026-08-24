set terminal pngcairo size 900,600
set output '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/validation_hpc_234/plots/trotter_S_pi_dq.png'
set title 'S_pi_dq vs dt^2 (correlated jackknife)'
set xlabel 'dt^2'
set ylabel 'S_pi_dq'
set key outside
f(x)=0.035286363700737046-0.0040440545683479989*x
plot '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/validation_hpc_234/plots/trotter_S_pi_dq.dat' using 1:2:3 with yerrorbars title 'QMC', f(x) title 'correlated GLS fit'
