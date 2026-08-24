set terminal pngcairo size 900,600
set output '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/validation_hpc_234/plots/blocking_L10_V2_existing.png'
set title 'Blocking errors: L10_V2_existing'
set xlabel 'block size'
set ylabel 'standard error'
set key outside
plot '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/validation_hpc_234/plots/blocking_L10_V2_existing.dat' using 1:2:3 with linespoints title 'sign', '' using 1:2:3 with linespoints title 'sign*S_pi', '' using 1:2:3 with linespoints title 'sign*S_pi_dq'
