set terminal pngcairo
set output 'sumTf.png'
set style data lines
plot 'globalSumTf_diff.dat' using 1:5 title "vol-weighted sum(Tf)"
