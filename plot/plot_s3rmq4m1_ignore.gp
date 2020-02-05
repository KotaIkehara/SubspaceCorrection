dt1 = "dt '.'"
dt2 = "dt '-'"
dt3 = "dt '_'"
dt4 = "dt '-.' lc 'royalblue'"
dt5 = "dt '_-'"
dt6  = "dt (1,4)"
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

plot   "../data/s3rmq4m1/ignore2/s3rmq4m1_ignore_zite=0_th=1.dat" using 1:2 with line title 'CG'
replot   "../data/s3rmq4m1/ignore2/s3rmq4m1_ignore_zite=1_th=1.dat" using 1:2 with line @dt1 title 'SC-CG({/Symbol a}=1, 2, 3, 4)'
replot   "../data/s3rmq4m1/ignore2/s3rmq4m1_ignore_zite=1_th=5.dat" using 1:2 with line @dt2 title 'SC-CG({/Symbol a}=5)'
replot   "../data/s3rmq4m1/ignore2/s3rmq4m1_ignore_zite=1_th=6.dat" using 1:2 with line @dt2 title 'SC-CG({/Symbol a}=6, 7, 8)'
replot   "../data/s3rmq4m1/ignore2/s3rmq4m1_ignore_zite=1_th=9.dat" using 1:2 with line @dt6 title 'SC-CG({/Symbol a}=9)'
set key outside