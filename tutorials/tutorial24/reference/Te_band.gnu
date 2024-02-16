set style data dots
set nokey
set xrange [0: 6.57944]
set yrange [ -0.72964 : 10.24675]
set arrow from  0.52728,  -0.72964 to  0.52728,  10.24675 nohead
set arrow from  1.60490,  -0.72964 to  1.60490,  10.24675 nohead
set arrow from  2.54472,  -0.72964 to  2.54472,  10.24675 nohead
set arrow from  3.35863,  -0.72964 to  3.35863,  10.24675 nohead
set arrow from  3.82853,  -0.72964 to  3.82853,  10.24675 nohead
set arrow from  4.35581,  -0.72964 to  4.35581,  10.24675 nohead
set arrow from  5.29563,  -0.72964 to  5.29563,  10.24675 nohead
set arrow from  6.10954,  -0.72964 to  6.10954,  10.24675 nohead
set xtics (" A "  0.00000," G "  0.52728," H "  1.60490," A "  2.54472," L "  3.35863," H "  3.82853," K "  4.35581," G "  5.29563," M "  6.10954," K "  6.57944)
 plot "Te_band.dat"
