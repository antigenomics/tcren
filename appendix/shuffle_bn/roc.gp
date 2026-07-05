# ROC curve — real vs shuffled Gaussian BN classifier (5-fold CV). Dark2 palette (Gnuplotting/ColorBrewer).
set terminal pdfcairo size 7.5cm,7.5cm font "Helvetica,11"
set output 'figures/roc.pdf'
set xlabel 'False positive rate'; set ylabel 'True positive rate'
set xrange [0:1]; set yrange [0:1]; set size square
set key bottom right box opaque; set grid lc rgb '#dddddd'
set style line 1 lc rgb '#1b9e77' lw 3
plot x w l lc rgb '#999999' dt 2 notitle, \
     'figures/roc.dat' u 1:2 w l ls 1 t 'BN classifier'
