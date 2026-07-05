# Posterior coefficient forest (standardised logistic weights, 94% credible interval), sorted.
set terminal pdfcairo size 11cm,15cm font "Helvetica,8" noenhanced
set output 'figures/forest.pdf'
set xlabel 'posterior coefficient (standardised, 94% CI)'
set grid xtics lc rgb '#dddddd'
set xzeroaxis lw 1 lc rgb '#cc3333'
set ytics font ",7"
set offsets 0,0,0.6,0.6
plot 'figures/forest.dat' u 2:1:3:4:ytic(5) w xerrorbars pt 7 ps 0.5 lc rgb '#1b9e77' notitle
