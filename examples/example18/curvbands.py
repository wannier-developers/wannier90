import pylab as pl
import numpy as np
import mpl_toolkits

ef=12.6128 # Fermi energy

## If False: save on outfile; if True: open window
outfile = 'Fe_curvbands.pdf'

######################################################################
##### There should be no need to modify anything below this line #####
######################################################################

x = np.loadtxt('fe_slice_x.dat')
dimx=x.size

y = np.loadtxt('fe_slice_y.dat')
dimy=y.size

z = np.loadtxt('fe_slice_curv.dat')
zz=z.reshape((dimx,dimy)).transpose()

bands=np.loadtxt('fe_slice_bands.dat')
numbands=bands.size/dimx/dimy
bands2=bands.reshape((dimx,dimy,numbands))

#title('Berry curvature of bcc Fe on the (010) plane')
for i in range(numbands):
    pl.contour(np.transpose(bands2[:,:,i]),[ef],colors='black')

pl.contourf(zz)
#imshow(zz,interpolation='bicubic',origin='lower')
pl.colorbar()

# Remove the axes
ax = pl.gca()
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

pl.savefig(outfile)
pl.show()
