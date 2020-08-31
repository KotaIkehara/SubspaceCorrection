dt1 = "lc '#a54675' lw 2.5"
dt2 = "lc '#296fbc' lw 2.5"
dt3 = "lc '#3d9435' lw 2.5"
dt4 = "lc '#cb360d' lw 2.5"
dt5 = "lc '#e1aa13' lw 2.5"
dt6  = "lc '#138bae' lw 2.5"
dt7  = "lc '#a5a5a5' lw 2.5"
dt8  = "dt (10,10)"
dt9  = "dt (20,10)"

set size square
set autoscale x
set autoscale y
set xlabel 'Iteration'
set ylabel 'Relative residual norm'
set logscale y
set format y " 10^{%T}"

plot   "../data2/rand/s3rmq4m1/ICCG/s3rmq4m1_iccg_zite=0_th=3.dat" using 1:2 with line @dt1 title  'ICCG'
replot   "../data/s3rmq4m1/IC2/s3rmq4m1_zite=1_th=5.dat" using 1:2 with line @dt3 title 'SC-ICCG({/Symbol a}=1, 2, 3, 4, 5)'
replot   "../data2/rand/s3rmq4m1/SCICCG/s3rmq4m1_sciccg_zite=1_th=6.dat" using 1:2 with line @dt3 title 'SC-ICCG({/Symbol a}=6)'
replot   "../data2/rand/s3rmq4m1/SCICCG/s3rmq4m1_sciccg_zite=1_th=7.dat" using 1:2 with line @dt3 title 'SC-ICCG({/Symbol a}=7)'
replot   "../data2/rand/s3rmq4m1/SCICCG/s3rmq4m1_sciccg_zite=1_th=8.dat" using 1:2 with line @dt3 title 'SC-ICCG({/Symbol a}=8)'
replot   "../data2/rand/s3rmq4m1/SCICCG/s3rmq4m1_sciccg_zite=1_th=9.dat" using 1:2 with line @dt3 title 'SC-ICCG({/Symbol a}=9)'
