set terminal wxt
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
set xyplane 0
set xzeroaxis
set yzeroaxis
set isosamples 6,6
set pm3d
set palette gray
unset colorbox
set style line 1 lt 1 lw 1
set style line 2 lt 2 lw 1
set style line 3 lt 2 pt 2 ps 3 lw 1
set style line 4 lt 3 lw 1 lc rgbcolor '#feb36c' pt 6 ps 3

splot "stencil36central.dat" using 1:2:3 with impulses ls 4 notitle, \
      "stencil36central.dat" using 1:2:3 with points ls 4 notitle, \
      "stencil36.dat" using 1:2:3 with impulses ls 2 notitle, \
      "stencil36.dat" using 1:2:3 with points ls 3 notitle
