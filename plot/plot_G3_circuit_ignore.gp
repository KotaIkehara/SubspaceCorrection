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

plot   "../data/G3_circuit/ignore2/G3_circuit_ignore_zite=0_th=1.dat" using 1:2 with line title 'CG'
replot   "../data/G3_circuit/ignore2/G3_circuit_ignore_zite=1_th=1.dat" using 1:2 with line @dt1 title 'SC-CG({/Symbol a}=1, 2, 3, 4)'

replot   "../data/G3_circuit/ignore2/G3_circuit_ignore_zite=1_th=5.dat" using 1:2 with line @dt2 title 'SC-CG({/Symbol a}=5)'