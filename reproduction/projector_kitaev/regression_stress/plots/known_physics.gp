set terminal pngcairo size 1000,650
set key outside
set xlabel 'V'
set ylabel 'R_cdw'
set output '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/regression_stress/plots/R_cdw_vs_V.png'
plot for [L in '10 18 26 34'] '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/regression_stress/plots/known_physics_R_vs_V.dat' using (strcol(1) eq L?$2:1/0):3:4 with yerrorlines title sprintf('L=%s',L)
set xlabel '(V-4)L (nu=1)'
set output '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/regression_stress/plots/R_cdw_collapse_nu1.png'
plot for [L in '10 18 26 34'] '/home/sunxr/PfQMC-main/reproduction/projector_kitaev/regression_stress/plots/known_physics_collapse_nu1.dat' using (strcol(1) eq L?$2:1/0):3:4 with yerrorlines title sprintf('L=%s',L)
