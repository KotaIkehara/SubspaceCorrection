dt1 = "lc '#a54675' lw 2.5"
dt2 = "lc '#296fbc' lw 2.5"
dt3 = "lc '#3d9435' lw 2.5"
dt4 = "lc '#cb360d' lw 2.5"
dt5 = "lc '#e1aa13' lw 2.5"
dt6  = "dt (1,4) lc '#138bae'"
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

plot     "../Hook_1498.mtx.iccg.zite1.theta3.0.thread1.dat" using 1:2 with line @dt1 title 'theta3'