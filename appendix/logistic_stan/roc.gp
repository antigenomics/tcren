# ROC — distribution-aware Bayesian logistic (real vs shuffled, 5-fold CV). Dark2 palette.
set terminal pdfcairo size 7.5cm,7.5cm font "Helvetica,11"
set output 'figures/roc.pdf'
set xlabel 'False positive rate'; set ylabel 'True positive rate'
set xrange [0:1]; set yrange [0:1]; set size square
set key bottom right box opaque; set grid lc rgb '#dddddd'
plot x w l lc rgb '#999999' dt 2 notitle, \
     'figures/roc.dat' u 1:2 w l lw 3 lc rgb '#1b9e77' t 'Bayesian logistic'
