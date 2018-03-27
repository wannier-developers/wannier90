import pylab as pl
import numpy as np
import matplotlib.mlab as ml
from collections import OrderedDict
 
points = np.loadtxt('iron_up-kslice-coord.dat')
points_x=points[:,0]
points_y=points[:,1]
num_pt=len(points)
 
area=    4.792898
 
x_coord=list(OrderedDict.fromkeys(points_x))
y_coord=list(OrderedDict.fromkeys(points_y))
dimx=len(x_coord)
dimy=len(y_coord)
 
# Energy level for isocontours (typically the Fermi level)
ef=   12.625600
 
bands_up=np.loadtxt('iron_up-kslice-bands.dat')
bands_dn=np.loadtxt('iron_dn-kslice-bands.dat')
numbands=bands_up.size//num_pt
bbands_up=bands_up.reshape((dimx,dimy,numbands))
bbands_dn=bands_dn.reshape((dimx,dimy,numbands))
for i in range(numbands):
  pl.contour(x_coord,y_coord,bbands_up[:,:,i],[ef],colors="red")
  pl.contour(x_coord,y_coord,bbands_dn[:,:,i],[ef],colors="blue")
 
# Remove the axes
ax = pl.gca()
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
 
pl.axes().set_aspect('equal')
 
outfile = 'iron_updn-fermi_lines.pdf'
 
 
pl.savefig(outfile)
pl.show()
