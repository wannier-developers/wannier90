import numpy as np
import matplotlib.pyplot as plt
import os

#Useful constants
j=complex(0,1)
eps0_SI=8.854187817e-12
hbar_SI=1.054571726e-34
elem_charge_SI=1.602176565e-19
hbar_eVs=6.582119e-16
bohr=0.52917721092
eV_au=3.674932379e-2
fac_SI2au=1e-10*eps0_SI*hbar_SI/(elem_charge_SI**2*hbar_eVs)*bohr/eV_au

cart=['x','y','z']

num_rows=np.zeros(9,dtype=int)
num_cols=np.zeros(9,dtype=int)

for a in range(3):
    for b in range(3):
        locals()['Na-eps_'+str(cart[a])+str(cart[b])]=np.loadtxt('Na-eps_'+str(cart[a])+str(cart[b])+'.dat')
        num_rows[a*3+b],num_cols[a*3+b]=locals()['Na-eps_'+str(cart[a])+str(cart[b])].shape

if all(i != num_rows[0] for i in num_rows):
    print('\n'+'Something wrong: the energy range do not coincide for some components of the linear KS tensors. Check it.'+'\n')
    exit_program()

energies=np.zeros((9,num_rows[0]))

for a in range(3):
    for b in range(3):
        for w in range(num_rows[0]):
            energies[a*3+b,w]=locals()['Na-eps_'+str(cart[a])+str(cart[b])][w,0]

for w in range(num_rows[0]):
    if all(ww != energies[0,w] for ww in energies[:,w]):
        print('\n'+'Something wrong: the values of the energy range do not coincide for some components of the linear KS tensors. Check it.'+'\n')
        exit_program()

alpha=np.zeros((3,3,num_rows[0]),dtype='complex_')
sigma=np.zeros((3,3,num_rows[0]),dtype='complex_')
eps=np.zeros((3,3,num_rows[0]),dtype='complex_')

for a in range(3):
    for b in range(3):
        for w in range(num_rows[0]):
            alpha[a,b,w]=locals()['Na-eps_'+str(cart[a])+str(cart[b])][w,1]+locals()['Na-eps_'+str(cart[a])+str(cart[b])][w,2]*j
            sigma[a,b,w]=locals()['Na-eps_'+str(cart[a])+str(cart[b])][w,3]+locals()['Na-eps_'+str(cart[a])+str(cart[b])][w,4]*j
        eps[a,b,:]=4.e0*np.pi*fac_SI2au*alpha[a,b,:]
    eps[a,a,:]=1.e0+eps[a,a,:]

epsi=np.zeros((3,3,num_rows[0]),dtype='complex_')

for w in range(num_rows[0]):
    epsi[:,:,w]=np.linalg.inv(eps[:,:,w])

format_string = '{:18.8E}'

for a in range(3):
    for b in range(3):
        with open('Na-eel_'+str(cart[a])+str(cart[b])+'.dat','w') as file:
            for w in range(num_rows[0]):
                numbers=np.array([energies[0,w],-epsi[a,b,w].imag])
                formatted_numbers = ' '.join([format_string.format(num) for num in numbers])
                file.write(formatted_numbers + '\n')
        file.close()
