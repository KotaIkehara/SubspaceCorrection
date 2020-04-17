dt1 = "lc '#a54675'"
dt2 = "lc '#296fbc'"
dt3 = "lc '#3d9435'"
dt4 = "lc '#cb360d'"
dt5 = "lc '#e1aa13'"
dt6  = "lc '#138bae'"
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
set yrange [0.000001:1000]

plot   "../data/s3rmq4m1/ignore2/s3rmq4m1_ignore_zite=0_th=1.dat" using 1:2 with line title 'CG'