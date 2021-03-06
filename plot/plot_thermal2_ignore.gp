dt1 = "lc '#a54675' lw 2.5"
dt2 = "lc '#296fbc' lw 2.5"
dt3 = "lc '#3d9435' lw 2.5"
dt4 = "lc '#cb360d' lw 2.5"
dt5 = "lc '#e1aa13' lw 2.5"
dt6  = "lc '#138bae' lw 2.5"
dt7  = "dt (5,5)"
dt8  = "dt (10,10)"
dt9  = "dt (20,10)"

set size square
set autoscale x
set autoscale y
set xlabel 'Iteration'
set ylabel 'Relative residual norm'
set logscale y
set format y " 10^{%T}"

plot   "../data2/rand/thermal2/ICCG/thermal2_iccg_zite=1_th=3.dat" using 1:2 with line @dt1 title 'ICCG'
replot   "../data2/rand/thermal2/SCICCG/thermal2_sciccg_zite=1_th=3.dat" using 1:2 with line @dt2 title 'SC-ICCG({/Symbol a}=3, 4)'
replot   "../data2/rand/thermal2/SCICCG/thermal2_sciccg_zite=1_th=5.dat" using 1:2 with line @dt3 title 'SC-ICCG({/Symbol a}=5)'
