# Precision-recall curve — real vs shuffled BN classifier (5-fold CV). Baseline = real fraction.
set terminal pdfcairo size 7.5cm,7.5cm font "Helvetica,11"
set output 'figures/pr.pdf'
base = real(system("cat figures/base.txt"))
set xlabel 'Recall'; set ylabel 'Precision'
set xrange [0:1]; set yrange [0:1]; set size square
set key bottom left box opaque; set grid lc rgb '#dddddd'
set style line 1 lc rgb '#d95f02' lw 3
plot base w l lc rgb '#999999' dt 2 t sprintf('baseline %.3f', base), \
     'figures/pr.dat' u 1:2 w l ls 1 t 'BN classifier'
