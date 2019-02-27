set style data dots
set nokey
set xrange [0: 6.57944]
set yrange [ -0.59933 :  9.98314]
set arrow from  0.52728,  -0.59933 to  0.52728,   9.98314 nohead
set arrow from  1.60490,  -0.59933 to  1.60490,   9.98314 nohead
set arrow from  2.54472,  -0.59933 to  2.54472,   9.98314 nohead
set arrow from  3.35863,  -0.59933 to  3.35863,   9.98314 nohead
set arrow from  3.82853,  -0.59933 to  3.82853,   9.98314 nohead
set arrow from  4.35581,  -0.59933 to  4.35581,   9.98314 nohead
set arrow from  5.29563,  -0.59933 to  5.29563,   9.98314 nohead
set arrow from  6.10954,  -0.59933 to  6.10954,   9.98314 nohead
set xtics (" A "  0.00000," G "  0.52728," H "  1.60490," A "  2.54472," L "  3.35863," H "  3.82853," K "  4.35581," G "  5.29563," M "  6.10954," K "  6.57944)
 plot "Te_band.dat"
