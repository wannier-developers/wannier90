#!/usr/bin/env python2
#
# gw2wannier90 interface
#
# This file is distributed as part of the Wannier90 code and
# under the terms of the GNU General Public License. See the
# file `LICENSE' in the root directory of the Wannier90
# distribution, or http://www.gnu.org/copyleft/gpl.txt
#
# The webpage of the Wannier90 code is www.wannier.org
#  
# The Wannier90 code is hosted on GitHub:
#
# https://github.com/wannier-developers/wannier90
# 
# Designed and tested with: Quantum Espresso and Yambo
# This interface should work with any G0W0 code 
# Originally written by Stepan Tsirkin
# Extended, developed and documented by Antimo Marrazzo
# 
# Last update February 19th, 2017 by Antimo Marrazzo (antimo.marrazzo@epfl.ch)
#
import numpy as np
import os,shutil
import datetime
from scipy.io import FortranFile
import sys
import glob


argv=sys.argv

if len(argv)<2:
    print "### gw2wannier90 interface ###"
    print "You need to provide the seedname"
    print "Usage: gw2wannier90.py seedname options"
    print "Options can be mmn, amn, spn, unk, uhu"
    print "Be careful with unformatted files, they are compiler-dependdent"
    print "A safer choice is to use (bigger) formatted files, in this case type:"
    print "spn_formatted, uiu_formatted, uhu_formatted" 
    print "If no options is specified, all the matrices and files are considered"    
    exit()
print '------------------------------\n'
print '##############################\n'
print '### gw2wannier90 interface ###\n'
print '##############################\n'

seedname=argv[1]   # for instance "silicon"
seednameGW=seedname+".gw"  #for instance "silicon.gw"
targets=[s.lower() for s in argv[2:]]   #options read from command line

#In case of formatted spn,uIu and uHu (mmn,amn,eig are formatted by default)
#NB: Formatted output is strongly reccommended! Fortran binaries are compilers dependent.
SPNformatted="spn_formatted" in targets   
UIUformatted="uiu_formatted" in targets
UHUformatted="uhu_formatted" in targets
UNKformatted="unk_formatted" in targets

if set(targets).intersection(set(["spn","uhu","mmn","amn","unk"])) :
    calcAMN="amn" in targets
    calcMMN="mmn" in targets
    calcUHU="uhu" in targets
    calcUIU="uiu" in targets
    calcSPN="spn" in targets
    calcUNK="unk" in targets
else:
    calcAMN=True
    calcMMN=True
    calcUHU=True
    calcUIU=True
    calcSPN=True
    calcUNK=True

if calcUHU : calcMMN=True

#Here we open a file to dump all the intermediate steps (mainly for debugging)
f_raw=open(seedname+'.gw2wannier90.raw','w')
#Opening seedname.nnkp file
f=open(seedname+".nnkp","r")
#It copies the seedname.win for GW, we should make this optional
#shutil.copy(seedname+".win",seednameGW+".win")
while True:
    s=f.readline()
    if "begin kpoints" in s: break
NKPT=int(f.readline())
print "Kpoints number:",NKPT
n1=np.array(NKPT,dtype=int)
IKP=[tuple(np.array(np.round(np.array(f.readline().split(),dtype=float)*n1),dtype=int)) for i in xrange (NKPT)]

while True:
    s=f.readline()
    if "begin nnkpts" in s: break
NNB=int(f.readline())
    
KPNB=np.array([[int(f.readline().split()[1])-1  for inb in xrange(NNB)] for ikpt in xrange(NKPT)])

while True:
    s=f.readline()
    if "begin exclude_bands" in s: break
exbands=np.array(f.readline().split(),dtype=int)
if len(exbands)>1 or exbands[0]!=0:
    #raise RuntimeError("exclude bands is not supported yet") # actually it is OK, see below
    print 'Exclude bands option is used: be careful to be consistent'
    'with the choice of bands for the GW QP corrections.'

corrections=np.loadtxt(seedname+".gw.unsorted.eig")
corrections={(int(l[1])-1,int(l[0])-1):l[2] for l in corrections}
print "G0W0 QP corrections read from ",seedname+".gw.unsorted.eig" 
#print corrections
eigenDFT=np.loadtxt(seedname+".eig")
nk=int(eigenDFT[:,1].max())
assert nk==NKPT

nbndDFT=int(eigenDFT[:,0].max())
eigenDFT=eigenDFT[:,2].reshape(NKPT,nbndDFT,order='C')
#print eigenDFT
f_raw.write('------------------------------\n')
f_raw.write('Writing DFT eigenvalues \n')
for line in eigenDFT:
    f_raw.write(str(line)+'\n')
f_raw.write('------------------------------ \n')

providedGW= [ ib for ib in xrange(nbndDFT) if all(  (ik,ib) in corrections.keys() for ik in xrange(NKPT)) ]
#print providedGW
f_raw.write('------------------------------\n')
f_raw.write('List of provided GW corrections (bands indeces)\n')
f_raw.write(str(providedGW))
f_raw.write('------------------------------\n')
NBND=len(providedGW)
print "Adding GW QP corrections to KS eigenvalues"
eigenDE= np.array([  [corrections[(ik,ib)]  for ib in providedGW ]   for ik in xrange(NKPT)])
eigenDFTGW= np.array([  [eigenDFT[ik,ib]+corrections[(ik,ib)]  for ib in providedGW ]   for ik in xrange(NKPT)])

f_raw.write('------------------------------\n')
f_raw.write('Writing GW eigenvalues unsorted (KS + QP correction) \n')
for line in eigenDFTGW:
    f_raw.write(str(line)+'\n')
f_raw.write('------------------------------\n')

print "Sorting"
bsort=np.array([np.argsort(eigenDFTGW[ik,:]) for ik in xrange(NKPT)])

f_raw.write('------------------------------\n')
f_raw.write('Writing sorting list\n')
for line in bsort:
    f_raw.write(str(line)+'\n')
f_raw.write('------------------------------\n')


eigenDE=np.array([eigenDE[ik][bsort[ik]] for ik in xrange(NKPT)])
eigenDFTGW=np.array([eigenDFTGW[ik][bsort[ik]] for ik in xrange(NKPT)])
BANDSORT=np.array([np.array(providedGW)[bsort[ik]] for ik in xrange(NKPT)])

f_raw.write('------------------------------\n')
f_raw.write('Writing sorted GW eigenvalues\n')
for line in eigenDFTGW:
    f_raw.write(str(line)+'\n')
f_raw.write('------------------------------\n')

print "GW eigenvalues sorted"
#print eigenDFT
print '------------------------------\n'
print "writing "+seednameGW+".eig"
feig_out=open(seednameGW+".eig","w")
for ik in xrange(NKPT):
    for ib in xrange(NBND):
	feig_out.write(" {0:4d} {1:4d} {2:17.12f}\n".format(ib+1,ik+1,eigenDFTGW[ik,ib]))
feig_out.close()
print seednameGW+".eig", ' written.'
('------------------------------\n')

if calcAMN:
  try:
    print "----------\n AMN module  \n---------\n"
    f_amn_out=open(seednameGW+".amn","w")
    f_amn_in=open(seedname+".amn","r")
    s=f_amn_in.readline().strip()
    f_amn_out.write("{0}, sorted by GW quasi-particle energies on {1} \n".format(s,datetime.datetime.now().isoformat()) )
    s=f_amn_in.readline()
    nb,nk,npr=np.array(s.split(),dtype=int)
    assert nk==NKPT
    assert nb==nbndDFT
    f_amn_out.write("  {0}   {1}    {2}   \n".format(NBND,nk,npr) )

    AMN=np.loadtxt(f_amn_in,dtype=float)[:,3:5]
    AMN=np.reshape(AMN[:,0]+AMN[:,1]*1j,(nb,npr,nk) ,order='F')
    for ik in xrange(nk):
	amn=AMN[BANDSORT[ik],:,ik]
	for ipr in xrange(npr):
	    for ib in xrange(NBND):
		f_amn_out.write(" {0:4d} {1:4d} {2:4d}  {3:16.12f}  {4:16.12f}\n".format(ib+1,ipr+1,ik+1,amn[ib,ipr].real,amn[ib,ipr].imag))
    f_amn_in.close()
    f_amn_out.close()
    print "----------\n AMN  - OK \n---------\n"
  except IOError as err:
    print "WARNING: {0}.amn not written : ".format(seednameGW) ,err


if calcMMN:
  try:
    print "----------\n MMN module  \n---------\n"

    f_mmn_out=open(os.path.join(seednameGW+".mmn"),"w")
    f_mmn_in=open(os.path.join(seedname+".mmn"),"r")

    s=f_mmn_in.readline().strip()
    f_mmn_out.write("{0}, sorted by GW quasi-particle energies on {1} \n".format(s,datetime.datetime.now().isoformat()) )
    s=f_mmn_in.readline()
    nb,nk,nnb=np.array(s.split(),dtype=int)
    assert nb==nbndDFT
    assert nk==NKPT
    f_mmn_out.write("    {0}   {1}    {2} \n".format(NBND,nk,nnb) )

    MMN=[]
    for ik in xrange(nk):
	MMN.append([])
	for ib in xrange(nnb):
	    s=f_mmn_in.readline()
	    f_mmn_out.write(s)
	    ik1,ik2=(int(i)-1 for i in s.split()[:2])
    	    assert(ik==ik1)
	    assert(KPNB[ik][ib]==ik2)
	    tmp=np.array([[f_mmn_in.readline().split() for m in xrange(nb)] for n in xrange(nb) ],dtype=str)
	    tmp=np.array(tmp[BANDSORT[ik2],:,:][:,BANDSORT[ik1],:],dtype=float)
	    tmp=(tmp[:,:,0]+1j*tmp[:,:,1]).T
	    MMN[ik].append(tmp)
	    for n in xrange(NBND): 
	        for m in xrange(NBND):
		    f_mmn_out.write( "  {0:16.12f}  {1:16.12f}\n".format(MMN[ik][ib][m,n].real,MMN[ik][ib][m,n].imag) )
    print "----------\n MMN OK  \n---------\n"
  except IOError as err:
    print "WARNING: {0}.mmn not written : ".format(seednameGW),err
    if calcUHU: 
	print "WARNING: {0}.uHu file also will not be written : ".format(seednameGW)
	calcUHU=False


def reorder_uXu(ext,formatted=False):
  try:
    print "----------\n  {0}  \n----------".format(ext)
    
    if formatted:
	f_uXu_in = open(seedname+"."+ext, 'r')
	f_uXu_out = open(seednameGW+"."+ext, 'w')
	header=f_uXu_in.readline()
	f_uXu_out.write(header)
	nbnd,NK,nnb=np.array(f_uXu_in.readline().split(),dtype=int)
	f_uXu_out.write("  ".join(str(x) for x in [NBND,NK,nnb])+"\n")
    else:
	f_uXu_in = FortranFile(seedname+"."+ext, 'r')
	f_uXu_out = FortranFile(seednameGW+"."+ext, 'w')
	header=f_uXu_in.read_record(dtype='c')
        f_uXu_out.write_record(header)
	nbnd,NK,nnb=np.array(f_uXu_in.read_record(dtype=np.int32))
	f_uXu_out.write_record(np.array([NBND,NK,nnb],dtype=np.int32))
    
    assert nbnd==nbndDFT
    print nbnd,NK,nnb
    
    if formatted:
	uXu=np.loadtxt(f_uXu_in).reshape(-1)
	start=0
	length=nbnd*nbnd

    for ik in xrange(NKPT):
	for ib2 in xrange(nnb):
	    for ib1 in xrange(nnb):
		if formatted:
		    A=uXu[start:start+length]
		    start+=length
		else:
		    A=f_uXu_in.read_record(dtype=np.complex)
		A=(A.reshape(nbnd,nbnd,order='F')[BANDSORT[KPNB[ik][ib2]],:][:,BANDSORT[KPNB[ik][ib1]]]+
			np.einsum('ln,lm,l->nm',MMN[ik][ib2].conj(),MMN[ik][ib1],eigenDE[ik]) ).reshape(-1,order='F')
		if formatted:
		    f_uXu_out.write("".join("{0:26.16e}  {1:26.16f}\n".format(x.real,x.imag) for x in A))
		else:
	    	    f_uXu_out.write_record(A)
    f_uXu_out.close()
    f_uXu_in.close()
    print "----------\n {0} OK  \n----------\n".format(ext)
  except IOError as err:
    print "WARNING: {0}.{1} not written : ".format(seednameGW,ext),err 

if calcUHU: reorder_uXu("uHu",UHUformatted)
if calcUIU: reorder_uXu("uIu",UIUformatted)


if calcSPN:
  try:
    print "----------\n SPN  \n---------\n"

    if SPNformatted:
	f_spn_in =  open(seedname+".spn", 'r')
	f_spn_out = open(seednameGW+".spn", 'w')
	header=f_spn_in.readline()
        f_spn_out.write(header)
        nbnd,NK=np.array(f_spn_in.readline().split(),dtype=np.int32)
	f_spn_out.write("  ".join(str(x) for x in (NBND,NKPT) ) )
    else:
	f_spn_in = FortranFile(seedname+".spn", 'r')
	f_spn_out = FortranFile(seednameGW+".spn", 'w')
	header=f_spn_in.read_record(dtype='c')
        f_spn_out.write_record(header)
        nbnd,NK=f_spn_in.read_record(dtype=np.int32)
	f_spn_out.write_record(np.array([NBND,NKPT],dtype=np.int32))

    print "".join(header)
    assert nbnd==nbndDFT


    indm,indn=np.tril_indices(nbnd)
    indmQP,indnQP=np.tril_indices(NBND)

    if SPNformatted:
	SPN=np.loadtxt(f_spn_in).reshape(-1)
	start=0
	length=(3*nbnd*(nbnd+1))/2

    for ik in xrange(NK):
	if SPNformatted:
	    A=SPN[start:start+length]
	    start+=length
	else:
	    A=np.zeros((3,nbnd,nbnd),dtype=np.complex)
	A[:,indn,indm]=f_spn_in.read_record(dtype=np.complex).reshape(3,nbnd*(nbnd+1)/2,order='F')
	A[:,indm,indn]=A[:,indn,indm].conj()
	check=np.einsum('ijj->',np.abs(A.imag))
	if check> 1e-10:
	    raise RuntimeError ( "REAL DIAG CHECK FAILED for spn: {0}".format(check) )
	A=A[:,:,BANDSORT[ik]][:,BANDSORT[ik],:][:,indnQP,indmQP].reshape((3*NBND*(NBND+1)/2),order='F')
	if formatted:
	    f_spn_out.write("".join("{0:26:16e} {1:26:16e}\n".format(x.real,x.imag) for x in A))
	else:
	    f_spn_out.write_record(A)

    f_spn_in.close()
    f_spn_out.close()
    print "----------\n SPN OK  \n---------\n"
  except IOError as err:
    print "WARNING: {0}.spn not written : ".format(seednameGW) ,err

if calcUNK:
    unkgwdir="UNK_GW"
    unkdftdir="UNK_DFT"
    files_list=[]
    for f_unk_name in  glob.glob("UNK*"):
        files_list.append(f_unk_name)

    try:
        os.mkdir(unkgwdir)
        os.mkdir(unkdftdir)
    except OSError:
        pass

    for f_unk_name in files_list:
      try:
	shutil.move('./'+f_unk_name,'./'+unkdftdir+'/')
        if UNKformatted:
	    f_unk_out=open(os.path.join(unkgwdir,f_unk_name),"w")
	    f_unk_in=open(os.path.join(unkdftdir,f_unk_name),"r")
	    nr1,nr2,nr3,ik,nbnd=np.array(f_unk_in.readline().split(),dtype=int)
	    NR=nr1*nr2*nr3
	    f_unk_out.write(" ".join(str(x) for x in (nr1,nr2,nr3,ik,NBND))+"\n")
	    f_unk_out.write("\n".join(
		    np.array([l.strip() for l in f_unk_in],dtype=str).reshape((nbnd,NR),order='C')[BANDSORT[ik-1],:].reshape(-1,order='C') ) )
	else:
	    f_unk_out=FortranFile(os.path.join(unkgwdir,f_unk_name),"w")
	    f_unk_in=FortranFile(f_unk_name,"r")
	    nr1,nr2,nr3,ik,nbnd=f_unk_in.read_record(dtype=np.int32)
	    f_unk_out.write_record(np.array([nr1,nr2,nr3,ik,NBND],dtype=np.int32))
	    unk=np.array([f_unk_in.read_record(dtype=np.complex) for ib in xrange(nbnd)] )[BANDSORT[ik-1],:]
	    for i in xrange(NBND):
		f_unk_out.write_record(unk[ib])
	f_unk_in.close()
        f_unk_out.close()
        shutil.move('./'+unkgwdir+'/'+f_unk_name,'./')
      except IOError as err:
        if err.errno==21:
    	    pass
    	else:
    	    raise err
    os.rmdir(unkgwdir)
    print 'UNK files are being reordered, old files coming from DFT are available in UNK_DFT'
f_raw.close()
