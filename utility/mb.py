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

#Some variables
cart=['x','y','z']
alpha_S = [1, 2, 3, 1, 1, 2]
beta_S = [1, 2, 3, 2, 3, 3]
valid_responses=['eps','shg','sc']
valid_kxc1_lrc_approx=['1','2','3']
valid_metal=['yes','no']

#Ask the user the seedname of files
seedname=input('\n'+'Which is the seedname of the calculations?'+'\n')

#Ask the user which optical response to work with
response=input('\n'+'Which optical response do you want to include many-body effects into? (Please choose one the following)'+'\n'+'eps: linear optical response'+'\n'+'shg: second-harmonic generation'+'\n'+'sc: shift current'+'\n'+'\n')

#Check if the answer is valid
if response not in valid_responses:
    print('\n'+'Invalid response. Please choose one of the following: eps/shg/sc'+'\n')
    quit()

#Initialize evaluation logicals
eval_shg=False
eval_sc=False
eval_eps=False

#Define the folder where files should be present
folder_path='./'

#Initialize the missing files array
missing_files=[]

#Check if there is any missing file for shg
if response == 'shg':

    #Check for the linear files
    for a in range(3):
        for b in range(3):
            file_name = seedname+'-eps_'+str(cart[a])+str(cart[b])+'.dat'
            file_path = os.path.join(folder_path, file_name)
            if not os.path.isfile(file_path):
                missing_files.append(file_name)

    #Check for the quadratic files
    for a in range(3):
        for bc in range(6):
            b = alpha_S[bc]-1
            c = beta_S[bc]-1
            file_name = seedname+'-shg_'+str(cart[a])+str(cart[b])+str(cart[c])+'.dat'
            file_path = os.path.join(folder_path, file_name)
            if not os.path.isfile(file_path):
                missing_files.append(file_name)

    #If any file is missing, say which and quit
    if missing_files:
        print('\n'+'The following files are missing:'+'\n')
        for file_name in missing_files:
            print(file_name)
        quit()
    else:
        print('\n'+'All needed files are present.'+'\n')

    #Set up the needed logicals
    eval_shg = True
    eval_eps = True

#Check if there is any missing file for sc
elif response == 'sc':

    #Check for the linear files
    for a in range(3):
        for b in range(3):
            file_name = seedname+'-eps_'+str(cart[a])+str(cart[b])+'.dat'
            file_path = os.path.join(folder_path, file_name)
            if not os.path.isfile(file_path):
                missing_files.append(file_name)

    #Check for the quadratic files
    for a in range(3):
        for bc in range(6):
            b = alpha_S[bc]-1
            c = beta_S[bc]-1
            file_name = seedname+'-sc_'+str(cart[a])+str(cart[b])+str(cart[c])+'.dat'
            file_path = os.path.join(folder_path, file_name)
            if not os.path.isfile(file_path):
                missing_files.append(file_name)

    #If any file is missing, say which and quit
    if missing_files:
        print('\n'+'The following files are missing:'+'\n')
        for file_name in missing_files:
            print(file_name)
        quit()
    else:
        print('\n'+'All needed files are present.'+'\n')

    #Set up the needed logicals
    eval_sc = True
    eval_eps = True

#Check if there is any missing file for eps
elif response == 'eps':

    #Check for the linear files
    for a in range(3):
        for b in range(3):
            file_name = seedname+'-eps_'+str(cart[a])+str(cart[b])+'.dat'
            file_path = os.path.join(folder_path, file_name)
            if not os.path.isfile(file_path):
                missing_files.append(file_name)

    #If any file is missing, say which and quit
    if missing_files:
        print('\n'+'The following files are missing:'+'\n')
        for file_name in missing_files:
            print(file_name)
        quit()
    else:
        print('\n'+'All needed files are present.'+'\n')

    #Set up the needed logicals
    eval_eps = True

#Read KS inputs. Check if all have consistent energy range and values within.
#1) Check if all the files have the same energy range, i.e. the same number of energy steps
if eval_eps:
    num_eps_rows=np.zeros(9,dtype=int)
    num_eps_cols=np.zeros(9,dtype=int)
    for a in range(3):
        for b in range(3):
            locals()['data_eps_'+str(cart[a])+str(cart[b])]=np.loadtxt(seedname+'-eps_'+str(cart[a])+str(cart[b])+'.dat')
            num_eps_rows[a*3+b],num_eps_cols[a*3+b]=locals()['data_eps_'+str(cart[a])+str(cart[b])].shape
    if not eval_shg and not eval_sc:
        num_rows=num_eps_rows

if eval_shg:
    num_shg_rows=np.zeros(18,dtype=int)
    num_shg_cols=np.zeros(18,dtype=int)
    for a in range(3):
        for bc in range(6):
            b=alpha_S[bc]-1
            c=beta_S[bc]-1
            locals()['data_shg_'+str(cart[a])+str(cart[b])+str(cart[c])]=np.loadtxt(seedname+'-shg_'+str(cart[a])+str(cart[b])+str(cart[c])+'.dat')
            num_shg_rows[a*6+bc],num_shg_cols[a*6+bc]=locals()['data_shg_'+str(cart[a])+str(cart[b])+str(cart[c])].shape
    num_rows=np.concatenate((num_eps_rows,num_shg_rows),axis=None)
        
if eval_sc:
    num_sc_rows=np.zeros(18,dtype=int)
    num_sc_cols=np.zeros(18,dtype=int)
    for a in range(3):
        for bc in range(6):
            b=alpha_S[bc]-1
            c=beta_S[bc]-1
            locals()['data_sc_'+str(cart[a])+str(cart[b])+str(cart[c])]=np.loadtxt(seedname+'-sc_'+str(cart[a])+str(cart[b])+str(cart[c])+'.dat')
            num_sc_rows[a*6+bc],num_sc_cols[a*6+bc]=locals()['data_sc_'+str(cart[a])+str(cart[b])+str(cart[c])].shape
    num_rows=np.concatenate((num_eps_rows,num_sc_rows),axis=None)

if all(i != num_rows[0] for i in num_rows):
    print('\n'+'Something wrong: the energy range do not coincide for some components of the linear KS tensors. Check it.'+'\n')
    exit_program()

#2) Check if all the values of the energy range coincide
if eval_eps:
    energies_eps=np.zeros((9,num_rows[0]))
    for a in range(3):
        for b in range(3):
            for w in range(num_rows[0]):
                energies_eps[a*3+b,w]=locals()['data_eps_'+str(cart[a])+str(cart[b])][w,0]
    if not eval_shg and not eval_sc:
        energies=np.zeros((9,num_rows[0]))
        for w in range(num_rows[0]):
             energies[:,w]=energies_eps[:,w]

if eval_shg:
    energies_shg=np.zeros((18,num_rows[0]))
    energies=np.zeros((27,num_rows[0]))
    for a in range(3):
        for bc in range(6):
            b=alpha_S[bc]-1
            c=beta_S[bc]-1
            for w in range(num_rows[0]):
                energies_shg[a*6+bc,w]=locals()['data_shg_'+str(cart[a])+str(cart[b])+str(cart[c])][w,0]
    for w in range(num_rows[0]):
        energies[:,w]=np.concatenate((energies_eps[:,w],energies_shg[:,w]),axis=None)

if eval_sc:
    energies_sc=np.zeros((18,num_rows[0]))
    energies=np.zeros((27,num_rows[0]))
    for a in range(3):
        for bc in range(6):
            b=alpha_S[bc]-1
            c=beta_S[bc]-1
            for w in range(num_rows[0]):
                energies_sc[a*6+bc,w]=locals()['data_sc_'+str(cart[a])+str(cart[b])+str(cart[c])][w,0]
    for w in range(num_rows[0]):
        energies[:,w]=np.concatenate((energies_eps[:,w],energies_sc[:,w]),axis=None)

for w in range(num_rows[0]):
    if all(ww != energies[0,w] for ww in energies[:,w]):
        print('\n'+'Something wrong: the values of the energy range do not coincide for some components of the linear KS tensors. Check it.'+'\n')
        exit_program()

print('\n'+'The energy range and the values within coincide for all the components of the linear KS tensors.'+'\n')

#Load the KS responses from postw90.x for working through the MB script
if eval_eps:
    alpha=np.zeros((3,3,num_rows[0]),dtype='complex_')
    sigma=np.zeros((3,3,num_rows[0]),dtype='complex_')
    eps=np.zeros((3,3,num_rows[0]),dtype='complex_')
    for a in range(3):
        for b in range(3):
            for w in range(num_rows[0]):
                alpha[a,b,w]=locals()['data_eps_'+str(cart[a])+str(cart[b])][w,1]+locals()['data_eps_'+str(cart[a])+str(cart[b])][w,2]*j
                sigma[a,b,w]=locals()['data_eps_'+str(cart[a])+str(cart[b])][w,3]+locals()['data_eps_'+str(cart[a])+str(cart[b])][w,4]*j
            eps[a,b,:]=4.e0*np.pi*fac_SI2au*alpha[a,b,:]
        eps[a,a,:]=1.e0+eps[a,a,:]

if eval_shg:
    shg_alpha=np.zeros((3,3,3,num_rows[0]),dtype='complex_')
    shg_sigma=np.zeros((3,3,3,num_rows[0]),dtype='complex_')
    for a in range(3):
        for bc in range(6):
            b=alpha_S[bc]-1
            c=beta_S[bc]-1
            for w in range(num_rows[0]):
                shg_alpha[a,b,c,w]=locals()['data_shg_'+str(cart[a])+str(cart[b])+str(cart[c])][w,1]+locals()['data_shg_'+str(cart[a])+str(cart[b])+str(cart[c])][w,2]*j
                shg_alpha[a,c,b,w]=shg_alpha[a,b,c,w]
                shg_sigma[a,b,c,w]=locals()['data_shg_'+str(cart[a])+str(cart[b])+str(cart[c])][w,3]+locals()['data_shg_'+str(cart[a])+str(cart[b])+str(cart[c])][w,4]*j
                shg_sigma[a,c,b,w]=shg_sigma[a,b,c,w]

if eval_sc:
    sc=np.zeros((3,3,3,num_rows[0]),dtype='complex_')
    for a in range(3):
        for bc in range(6):
            b=alpha_S[bc]-1
            c=beta_S[bc]-1
            for w in range(num_rows[0]):
                sc[a,b,c,w]=locals()['data_sc_'+str(cart[a])+str(cart[b])+str(cart[c])][w,1]
                sc[a,c,b,w]=sc[a,b,c,w]

#Inform the user about the tensorial LRC xc kernel
print('\n'+'In the long-wavelength and optical limit, the tensorial exchange-correlation kernel takes the form of a long-range contribution (LRC) that we simplify to an attractive and static Coulomb-like potential.'+'\n')

#Define the LRC tensorial xc kernel
kxc1_LRC=np.zeros((3,3,num_rows[0]),dtype='complex_')

#Ask the user which optical response to work with
kxc1_lrc_approx=input('\n'+'The coefficients governing this potential can be defined by:'+'\n'+'1) A non-selfconsistent bootstrap calculation, also known as RPA-bootstrap (RBO)'+'\n'+'2) A self-consistent bootstrap calculation (BO)'+'\n'+'3) Empirical values'+'\n'+'\n'+'Which procedure do you want to use for defining the tensorial LRC xc kernel? (Please, give only the number)'+'\n'+'\n')

#Check if the answer is valid
if kxc1_lrc_approx not in valid_kxc1_lrc_approx:
    print('\n'+'Invalid response'+'\n')
    quit()

#Initialize evaluation logicals
eval_rbo=False
eval_bo=False
eval_emp=False

if kxc1_lrc_approx == '1':
    eval_rbo=True
elif kxc1_lrc_approx == '2':
    eval_bo=True
elif kxc1_lrc_approx == '3':
    eval_emp=True

#Bootstrap procedures (both RBO and BO)
if eval_rbo or eval_bo:

    metal=input('\n'+'Is the material metallic? yes/no'+'\n')
    if metal not in valid_metal:
        print('\n'+'Invalid response'+'\n')
        quit()

    if metal=='yes':
        print('\n'+'Since the material is metallic, the bootstrap approximations do not work. You should choose the option 3.'+'\n')
        quit()

    kxc1_BO=np.zeros((3,3,1000),dtype='complex_')
    lkxc1_BO=np.zeros((9,1000),dtype='complex_')

    epsi_BO=np.zeros((3,3,num_rows[0]),dtype='complex_')
    epsim_BO=np.zeros((3,3,num_rows[0]),dtype='complex_')
    fac_BO=-np.linalg.inv(fac_SI2au*alpha[:,:,0])

    for w in range(num_rows[0]):
        epsi_BO[:,:,w]=np.linalg.inv(eps[:,:,w])

    #selfconsistent loop with index i
    for i in range(1000):
        for w in range(num_rows[0]):
            for a in range(3):
                for b in range(3):
                    epsim_BO[a,b,w]=-4.e0*np.pi*sum(fac_SI2au*alpha[a,:,w]*epsi_BO[:,b,w])
                epsim_BO[a,a,w]=1.e0+epsim_BO[a,a,w]
        for a in range(3):
            for b in range(3):
                kxc1_BO[a,b,i]=sum(fac_BO[a,:]*epsim_BO[:,b,0])
                lkxc1_BO[a*3+b,i]=sum(fac_BO[a,:]*epsim_BO[:,b,0])
        nc=i
        #If RBO and i=0, the loop stops here: non-selfconsistent
        if eval_rbo and i==0:
            break
        #If BO and i>0, continue looping until convergence through self-consistency
        elif eval_bo and i>0:
            if all(diff < 1.e-8 for diff in abs(lkxc1_BO[:,i]-lkxc1_BO[:,i-1])):
                break
        for w in range(num_rows[0]):
            for a in range(3):
                for b in range(3):
                    epsi_BO[a,b,w]=fac_SI2au*(4.e0*np.pi*alpha[a,b,w]+sum(kxc1_BO[a,:,i]*alpha[:,b,w]))
                epsi_BO[a,a,w]=1.e0+epsi_BO[a,a,w]
            epsi_BO[:,:,w]=np.linalg.inv(epsi_BO[:,:,w])

    for w in range(num_rows[0]):
        kxc1_LRC[:,:,w]=kxc1_BO[:,:,nc]

#Empircal procedure
elif eval_emp:

    kxc1_emp=np.zeros((3,3),dtype='complex_')
    lkxc1_emp=np.zeros(9,dtype='complex_')

    lkxc1_emp_str=input('\n'+'Enter the empirical values you want for the tensorial LRC xc kernel K^{ab}_{\text{x}c,1} with numbers separated by commas, as: xx, xy, xz, yx, yy, yz, zx, zy, zz'+'\n'+'\n')

    # Split the input string into individual values, and convert them to integers or floats
    try:
        lkxc1_emp=[int(value) if value.isdigit() else float(value) for value in lkxc1_emp_str.split(',')]
    except ValueError:
        print('Invalid input. Please enter numbers separated by commas.')
        quit()

    for a in range(3):
        for b in range (3):
            kxc1_emp[a,b]=lkxc1_emp[a*3+b]

    for w in range(num_rows[0]):
        kxc1_LRC[:,:,w]=kxc1_emp[:,:]

#Print the tensorial LRC xc kernel
if eval_rbo:
    print('\n'+'RPA-bootstrap approximation for the LRC tensorial xc kernel achieved, with values'+'\n')
elif eval_bo:
    print('\n'+'Bootstrap approximation for the LRC tensorial xc kernel achieved in',nc,'cycles, with values:'+'\n')
elif eval_emp:
    print('\n'+'The empirical values given for the LRC tensorial xc kernel are:'+'\n')

print(np.array_str(kxc1_LRC[:,:,0], precision=3, suppress_small=True))

#Calculate microscopic responses
if eval_eps:
    eps_micro=np.zeros((3,3,num_rows[0]),dtype='complex_')
    epsi_micro=np.zeros((3,3,num_rows[0]),dtype='complex_')
    sigma_micro=np.zeros((3,3,num_rows[0]),dtype='complex_')
    alpha_micro=np.zeros((3,3,num_rows[0]),dtype='complex_')
    for w in range(num_rows[0]):
        for a in range(3):  
            for b in range(3):
                eps_micro[a,b,w]=fac_SI2au*(4.e0*np.pi*alpha[a,b,w]+sum(kxc1_LRC[a,:,w]*alpha[:,b,w]))
            eps_micro[a,a,w]=1.e0+eps_micro[a,a,w]
        epsi_micro[:,:,w]=np.linalg.inv(eps_micro[:,:,w])
        for a in range(3):
            for b in range(3):
                alpha_micro[a,b,w]=sum(alpha[a,:,w]*epsi_micro[:,b,w])
                sigma_micro[a,b,w]=sum(sigma[a,:,w]*epsi_micro[:,b,w])

if eval_shg:
    shg_sigma_micro=np.zeros((3,3,3,(num_rows[0]-1)//2+1),dtype='complex_')
    shg_alpha_micro=np.zeros((3,3,3,(num_rows[0]-1)//2+1),dtype='complex_')
    for w in range((num_rows[0]-1)//2+1):
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    for d in range(3):
                        for e in range(3):
                            for f in range(3):
                                shg_alpha_micro[a,b,c,w]=shg_alpha_micro[a,b,c,w]+epsi_micro[a,d,2*w]*shg_alpha[d,e,f,w]*epsi_micro[e,b,w]*epsi_micro[f,c,w]
                                shg_sigma_micro[a,b,c,w]=shg_sigma_micro[a,b,c,w]+epsi_micro[a,d,2*w]*shg_sigma[d,e,f,w]*epsi_micro[e,b,w]*epsi_micro[f,c,w]

if eval_sc:
    sc_micro=np.zeros((3,3,3,num_rows[0]),dtype='complex_')
    for w in range(num_rows[0]):
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    for d in range(3):
                        for e in range(3):
                            for f in range(3):
                                sc_micro[a,b,c,w]=sc_micro[a,b,c,w]+epsi_micro[a,d,0]*sc[d,e,f,w]*epsi_micro[e,b,w]*np.conjugate(epsi_micro[f,c,w])

#Calculate macroscopic responses
if eval_eps:
    eps_macro=np.zeros((3,3,num_rows[0]),dtype='complex_')
    epsi_macro=np.zeros((3,3,num_rows[0]),dtype='complex_')
    sigma_macro=np.zeros((3,3,num_rows[0]),dtype='complex_')
    alpha_macro=np.zeros((3,3,num_rows[0]),dtype='complex_')
    for w in range(num_rows[0]):
        for a in range(3):
            for b in range(3):
                epsi_macro[a,b,w]=-fac_SI2au*4.e0*np.pi*alpha_micro[a,b,w]
            epsi_macro[a,a,w]=1.e0+epsi_macro[a,a,w]
        eps_macro[:,:,w]=np.linalg.inv(epsi_macro[:,:,w])
        for a in range(3):
            for b in range(3):
                alpha_macro[a,b,w]=sum(alpha_micro[a,:,w]*eps_macro[:,b,w])
                sigma_macro[a,b,w]=sum(sigma_micro[a,:,w]*eps_macro[:,b,w])

if eval_shg:
    shg_sigma_macro=np.zeros((3,3,3,(num_rows[0]-1)//2+1),dtype='complex_')
    shg_alpha_macro=np.zeros((3,3,3,(num_rows[0]-1)//2+1),dtype='complex_')
    for w in range((num_rows[0]-1)//2+1):
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    for d in range(3):
                        for e in range(3):
                            for f in range(3):
                                shg_alpha_macro[a,b,c,w]=shg_alpha_macro[a,b,c,w]+eps_macro[a,d,2*w]*shg_alpha_micro[d,e,f,w]*eps_macro[e,b,w]*eps_macro[f,c,w]
                                shg_sigma_macro[a,b,c,w]=shg_sigma_macro[a,b,c,w]+eps_macro[a,d,2*w]*shg_sigma_micro[d,e,f,w]*eps_macro[e,b,w]*eps_macro[f,c,w]

if eval_sc:
    sc_macro=np.zeros((3,3,3,num_rows[0]),dtype='complex_')
    for w in range(num_rows[0]):
        for a in range(3):
            for b in range(3):
                for c in range(3):
                    for d in range(3):
                        for e in range(3):
                            for f in range(3):
                                sc_macro[a,b,c,w]=sc_macro[a,b,c,w]+eps_macro[a,d,0]*sc_micro[d,e,f,w]*eps_macro[e,b,w]*np.conjugate(eps_macro[f,c,w])

#Write macroscopic (including MB effects) responses
format_string = '{:18.8E}'

if eval_eps:
    for a in range(3):
        for b in range(3):
            with open(seedname+'-eps-mb_'+str(cart[a])+str(cart[b])+'.dat','w') as file:
                for w in range(num_rows[0]):
                    numbers=np.array([energies[0,w],alpha_macro[a,b,w].real,alpha_macro[a,b,w].imag,sigma_macro[a,b,w].real,sigma_macro[a,b,w].imag,eps_macro[a,b,w].real,eps_macro[a,b,w].imag])
                    formatted_numbers = ' '.join([format_string.format(num) for num in numbers])
                    file.write(formatted_numbers + '\n')
            file.close()

if eval_shg:
    for a in range(3):
        for bc in range(6):
            b=alpha_S[bc]-1
            c=beta_S[bc]-1
            with open(seedname+'-shg-mb_'+str(cart[a])+str(cart[b])+str(cart[c])+'.dat','w') as file:
                for w in range((num_rows[0]-1)//2+1):
                    numbers=np.array([energies[0,w],shg_alpha_macro[a,b,c,w].real,shg_alpha_macro[a,b,c,w].imag,shg_sigma_macro[a,b,c,w].real,shg_sigma_macro[a,b,c,w].imag])
                    formatted_numbers = ' '.join([format_string.format(num) for num in numbers])
                    file.write(formatted_numbers + '\n')
            file.close()

if eval_sc:
    for a in range(3):
        for bc in range(6):
            b=alpha_S[bc]-1
            c=beta_S[bc]-1
            with open(seedname+'-sc-mb_'+str(cart[a])+str(cart[b])+str(cart[c])+'.dat','w') as file:
                for w in range(num_rows[0]):
                    numbers=np.array([energies[0,w],sc_macro[a,b,c,w].real,sc_macro[a,b,c,w].imag])
                    formatted_numbers = ' '.join([format_string.format(num) for num in numbers])
                    file.write(formatted_numbers + '\n')
            file.close()
