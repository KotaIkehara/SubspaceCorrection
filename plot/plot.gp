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
set xlabel 'iteration'
set ylabel 'residual'
set logscale y
plot   "data/s3rmq4m1/data_0.dat" using 1:2 with line title 'zite=0'
replot   "data/s3rmq4m1/IC/s3rmq4m1_zite=1_th=1.dat" using 1:2 with line @dt1 title 'm=10, θ=10^{-1}, 10^{-2}, 10^{-3}, 10^{-4}'
replot   "data/s3rmq4m1/IC/s3rmq4m1_zite=1_th=5.dat" using 1:2 with line @dt2 title 'm=8, θ=10^{-5}'
replot   "data/s3rmq4m1/IC/s3rmq4m1_zite=1_th=6.dat" using 1:2 with line @dt3 title 'm=7, θ=10^{-6}'
replot   "data/s3rmq4m1/IC/s3rmq4m1_zite=1_th=7.dat" using 1:2 with line @dt4 title 'm=6, θ=10^{-7}'
replot   "data/s3rmq4m1/IC/s3rmq4m1_zite=1_th=8.dat" using 1:2 with line @dt5 title 'm=5, θ=10^{-8}'
replot   "data/s3rmq4m1/IC/s3rmq4m1_zite=1_th=9.dat" using 1:2 with line @dt6 title 'm=4, θ=10^{-9}'
