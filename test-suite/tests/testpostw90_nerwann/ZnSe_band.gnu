set style data dots
set nokey
set xrange [0: 3.40668]
set yrange [ -0.05737 : 12.58560]
set arrow from  0.78385,  -0.05737 to  0.78385,  12.58560 nohead
set arrow from  1.74387,  -0.05737 to  1.74387,  12.58560 nohead
set arrow from  2.85241,  -0.05737 to  2.85241,  12.58560 nohead
set xtics ("W"  0.00000,"L"  0.78385,"G"  1.74387,"X"  2.85241,"W"  3.40668)
 plot "ZnSe_band.dat"
