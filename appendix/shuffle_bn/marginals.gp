# Marginals of the 6 most class-separating features: real (green) vs shuffled (orange).
set terminal pdfcairo size 18cm,11cm font "Helvetica,10"
set output 'figures/marginals.pdf'
set multiplot layout 2,3
set style fill transparent solid 0.45 noborder
set key top right font ",7"
do for [i=0:5] {
  set title sprintf("feature %d", i+1) font ",9"
  unset ylabel; set ytics format ""
  plot 'figures/marginals.dat' index i u 1:2 w filledcurves y1=0 lc rgb '#1b9e77' t 'real', \
       '' index i u 1:3 w filledcurves y1=0 lc rgb '#d95f02' t 'shuffled'
}
unset multiplot
