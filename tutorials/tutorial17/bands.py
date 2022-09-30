import pylab as pl
import numpy as np
import mpl_toolkits

ef=12.6128 # Fermi energy

x = np.loadtxt('fe_slice_x.dat')
dimx=x.size
y = np.loadtxt('fe_slice_y.dat')
dimy=y.size

bands=np.loadtxt('fe_slice_bands.dat')
numbands=bands.size/dimx/dimy
bands2=bands.reshape((dimx,dimy,numbands))

#title('Fermi contours of bcc Fe on the (010) plane')
for i in range(numbands):
    pl.contour(np.transpose(bands2[:,:,i]),[ef],colors='black')


# Remove the axes
ax = pl.gca()
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

pl.axes().set_aspect('equal')

#pl.savefig('Fe_fermicontours.png')
pl.savefig('Fe_fermicontours.pdf')
pl.show()
