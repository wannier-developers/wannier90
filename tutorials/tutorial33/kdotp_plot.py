import numpy as N
import sys   as SYS
import os    as OS

import matplotlib
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties

from matplotlib import rc
rc('text', usetex=True)

font       = FontProperties()
font.set_size(20)

#------------------------------------------------------
#       Extract W90 bands data
#------------------------------------------------------

W_band_dir   = './'
W_band_file  = 'bc2n_band.dat'

eF=-1.6103

data  = N.loadtxt(W_band_dir+W_band_file)

list_xW = data[:,0] 
bandsW       = data[:,1:].transpose()-eF 


#----------------------------------------------------------
# plot
#----------------------------------------------------------

fig = plt.figure(figsize=(10,7))

axplot  = fig.add_subplot(111)

for band in bandsW:
  axplot.plot(list_xW,band,linestyle='',label='ab initio (W90)', color='k',marker='o',markersize=2)



######################################
######################################
#######  extract k.p coefficients ####

seed = 'bc2n'

filename_0 =  seed + '-kdotp_0.dat'
filename_1 =  seed + '-kdotp_1.dat'
filename_2 =  seed + '-kdotp_2.dat'

# arrays containing pauli matrix coefficients proportional to 
# kx, ky, kx**2, ky**2 and kx*ky, 5 entries in total
a1 = [0] * 5
a2 = [0] * 5
a3 = [0] * 6 # contains extra term independent of k
a0 = [0] * 5

# assign spatial directions
# indexes: 0=x, 1=y, 2=z
index_x = 2
index_y = 0 


#### 0th order ####
data	= N.loadtxt(filename_0)
re = data[:,0]
im = data[:,1]

# Az
a3[0] = 0.5*N.real(re[0]-re[3])
offset = re[0] - a3[0]


#### 1st order ####
data	= N.loadtxt(filename_1)
re = data[:,0]
im = data[:,1]
full = re + im*1.j

# Az
# index_i*4 chooses the matrix chunk corresponding to the proper direction (the 4 factor comes from 2x2 bands), 
# +1 gets the band-off-diagonal entry
a3[1] =  0.5*(full[index_x*4]-full[index_x*4+3]).real
a3[2] =  0.5*(full[index_y*4]-full[index_y*4+3]).real

# Ax
# index_i*4 chooses the matrix chunk corresponding to the proper direction (the 4 factor comes from 2x2 bands), 
# +1 gets the band-off-diagonal entry
a1[0] =  (full[index_x*4+1]).real
a1[1] = -(full[index_y*4+1]).imag

# Ay
# index_i*4 chooses the matrix chunk corresponding to the proper direction (the 4 factor comes from 2x2 bands), 
# +1 gets the band-off-diagonal entry
a2[0] = -(full[index_x*4+1]).imag
a2[1] =  (full[index_y*4+1]).real

# A0
# index_i*4 chooses the matrix chunk corresponding to the proper direction (the 4 factor comes from 2x2 bands), 
# +1 gets the band-off-diagonal entry
a0[0] =  0.5*(full[index_x*4]+full[index_x*4+3]).real
a0[1] =  0.5*(full[index_y*4]+full[index_y*4+3]).real



#### 2nd order ####
data	= N.loadtxt(filename_2)
re = data[:,0]
im = data[:,1]
full = re + im*1.j

# a1
# first index jumps 12=3*4 (3 for spatial indexes, 4=2x2 for bands)
# second index only jumps 4=2x2 for bands
# +1 gets the band-off-diagonal term
a1[2] =  full[index_x*(12+4)+1].real 
a1[3] = -full[index_y*(12+4)+1].imag 
a1[4] =  full[index_x*12+index_y*4+1].real + full[index_y*12+index_x*4+1].real 

# a2
# first index jumps 12=3*4 (3 for spatial indexes, 4=2x2 for bands)
# second index only jumps 4=2x2 for bands
# +1 gets the band-off-diagonal term
a2[2] = -full[index_x*(12+4)+1].imag 
a2[3] =  full[index_y*(12+4)+1].real 
a2[4] = -full[index_x*12+index_y*4+1].imag - full[index_y*12+index_x*4+1].imag 

# a3
# first index jumps 12=3*4 (3 for spatial indexes, 4=2x2 for bands)
# second index only jumps 4=2x2 for bands
# +0 and +3 get the 11 and 22 band-diagonal terms 
a3[3] = N.real(0.5*(full[index_x*(12+4)] - full[index_x*(12+4)+3]))
a3[4] = N.real(0.5*(full[index_y*(12+4)] - full[index_y*(12+4)+3]))
a3[5] = 0.5*N.real( (full[index_x*12+index_y*4]+full[index_y*12+index_x*4] ) - (full[index_x*12+index_y*4+3] + full[index_y*12+index_x*4+3] ) )

# a0
a0[2] = N.real(0.5*(full[index_x*(12+4)] + full[index_x*(12+4)+3]))
a0[3] = N.real(0.5*(full[index_y*(12+4)] + full[index_y*(12+4)+3]))
a0[4] = 0.5*N.real( (full[index_x*12+index_y*4]+full[index_y*12+index_x*4] ) + (full[index_x*12+index_y*4+3] + full[index_y*12+index_x*4+3] ) )


# set up k.p band dispersion
# energy is in eV
# linear coefficients are in units of eV*Ang 
# quadratic coefficients are in units of eV*Ang**2 
a0_angs = 1.0
bxfac = -0.5*0.231321    
byfac = 0.0
bfac = N.sqrt( (bxfac)**(2) + (byfac)**(2) ) 
cos_xangle = bxfac/bfac
sin_yangle = byfac/bfac
kfac = bfac*2*N.pi/a0_angs

# define quantities to be plotted
kk_list = []
band_1  = []
band_2  = []
band_1lin  = []
band_2lin  = []

# k along x 
num_k = 200
for ii in range(0,num_k+1):  

  kk  = kfac*float(ii)/float(num_k)    
  kkx = cos_xangle*kk 
  kky = sin_yangle*kk 

  # coefficients with quadratic dispersion
  a1_coeff = (a1[0])*(kkx) + (a1[1])*(kky) + (a1[2])*(kkx**(2)) + (a1[3])*(kky**(2)) + (a1[4])*(kky*kkx) 
  a2_coeff = (a2[0])*(kkx) + (a2[1])*(kky) + (a2[2])*(kkx**(2)) + (a2[3])*(kky**(2)) + (a2[4])*(kky*kkx) 
  a3_coeff = a3[0] + (a3[1])*(kkx) + (a3[2])*(kky) + (a3[3])*(kkx**(2)) + (a3[4])*(kky**(2)) + (a3[5])*(kky*kkx)
  a0_coeff = a0[0]*(kkx) + a0[1]*(kky) + a0[2]*(kkx)**(2) + a0[3]*(kky)**(2) + (a0[4])*(kky*kkx)

  eps = N.sqrt((N.abs(a1_coeff-a2_coeff*(1.j)))**(2) + (N.abs(a3_coeff))**(2)) 
  
  kk_list.append(kk)
  band_1.append(a0_coeff+eps+offset-eF)
  band_2.append(a0_coeff-eps+offset-eF)

  # coefficients with linear dispersion
  a1_coeff = (a1[0])*(kkx) + (a1[1])*(kky) 
  a2_coeff = (a2[0])*(kkx) + (a2[1])*(kky) 
  a3_coeff = a3[0] + (a3[1])*(kkx) + (a3[2])*(kky) 
  a0_coeff = a0[0]*(kkx) + a0[1]*(kky) 

  eps = N.sqrt((N.abs(a1_coeff-a2_coeff*(1.j)))**(2) + (N.abs(a3_coeff))**(2)) 
  
  band_1lin.append(a0_coeff+eps+offset-eF)
  band_2lin.append(a0_coeff-eps+offset-eF)



axplot.plot(kk_list,band_1lin,color='b',linestyle='--',label='linear $k\cdot p$')
axplot.plot(kk_list,band_2lin,color='b',linestyle='--')

axplot.plot(kk_list,band_1,color='r',linestyle='-',label='quadratic $k\cdot p$')
axplot.plot(kk_list,band_2,color='r',linestyle='-')



######################################
######################################
## set up figure

axplot.set_xlim([list_xW[0],list_xW[-1]])

xtick_positions = [0.0,list_xW[len(list_xW)-1]]
xtick_labels = ['$\mathrm{S}$','$\mathrm{X}$']

axplot.set_xticks(xtick_positions)
axplot.set_xticklabels(xtick_labels,fontproperties=font)

for y_tick in axplot.yaxis.get_major_ticks():
    y_tick.label1.set_fontsize(20)

#       Set the yticks properties
axplot.set_ylabel('$E-E_{F}$ $\mathrm{(eV)}$',fontproperties=font)


axplot.axhline(0, color='k',linestyle='-',lw=0.25)
axplot.legend(loc=0,fancybox=True,shadow=True,prop=font)

fig.tight_layout()

plt.savefig('kdotp_bands_SX.pdf')
plt.show()
