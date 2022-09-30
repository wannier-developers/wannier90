#!/usr/bin/env gnuplot
set terminal pdf enhanced color dashed lw 1 size 6in,9in
set output "graphene_bandsdiff.pdf"
set size ratio 1.5
fermi = -2.3043
# set style data dots
set key
set xrange [0: 4.02877]
set yrange [(-22.91167 - fermi) : (23.53815 - fermi)]
set arrow from  1.47463, (-22.91167 - fermi) to  1.47463,  (23.53815 - fermi) nohead
set arrow from  2.32601, (-22.91167 - fermi) to  2.32601,  (23.53815 - fermi) nohead
set xtics ("G"  0.00000,"M"  1.47463,"K"  2.32601,"G"  4.02877)
# scale QE x-axis to be consistent with w90
plot "graphene.bands.dat.gnu" u ($1*2.554):($2 - fermi) w l title "QE", \
     "graphene_band.dat" u 1:($2 - fermi) w l title "W90"
