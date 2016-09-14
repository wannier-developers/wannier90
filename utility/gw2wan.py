#!/usr/bin/env python2


# The interfac between Quantum Espresso, (oeriginally for YAMBO) and wannier90
# Written by Stepan Tsirkin
#   12-16/09/2016,  Wannier90  workshop,  Donostia-Sansebastian.
# this version is not tested yet !


import netCDF4
import numpy as np
import os,shutil
from netCDF4 import Dataset
import datetime
from scipy.io import FortranFile
import sys


argv=sys.argv

if len(argv)<2:
    print "you need to provide the seedname"
    exit()

seedname=argv[1]
seednameGW=seedname+"+GW"
targets=[s.lower() for s in argv[2:]]

SPNformatted="spn_formatted" in targets
UIUformatted="uiu_formatted" in targets
UHUformatted="uhu_formatted" in targets

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

f=open(seedname+".nnkp","r")
shutil.copy(seedname+".nnkp",seednameGW+".nnkp")
while True:
    s=f.readline()
    if "begin kpoints" in s: break
NKPT=int(f.readline())
print "NK=",NKPT
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
    raise RuntimeError("exclude bands is not supported yet")

corrections=np.loadtxt(seedname+".gw.eig")
corrections={(int(l[1])-1,int(l[0])-1):l[2] for l in corrections}
print "eig-read"
#print corrections
eigenDFT=np.loadtxt(seedname+".eig")
nk=int(eigenDFT[:,1].max())
assert nk==NKPT

nbndDFT=int(eigenDFT[:,0].max())
eigenDFT=eigenDFT[:,2].reshape(NKPT,nbndDFT,order='C')
print eigenDFT
print "check"
providedGW= [ ib for ib in xrange(nbndDFT) if all(  (ik,ib) in corrections.keys() for ik in xrange(NKPT)) ]
print providedGW
NBND=len(providedGW)
print "adding"
eigenDE= np.array([  [corrections[(ik,ib)]  for ib in providedGW ]   for ik in xrange(NKPT)])
eigenDFTGW= np.array([  [eigenDFT[ik,ib]+corrections[(ik,ib)]  for ib in providedGW ]   for ik in xrange(NKPT)])
#print eigenDFTGW
print "eigenDFTGW -unsorted"
print eigenDFTGW

print "sorting"
bsort=np.array([np.argsort(eigenDFTGW[ik,:]) for ik in xrange(NKPT)])
print "bsort="
print bsort
eigenDE=np.array([eigenDE[ik][bsort[ik]] for ik in xrange(NKPT)])
eigenDFTGW=np.array([eigenDFTGW[ik][bsort[ik]] for ik in xrange(NKPT)])
BANDSORT=np.array([np.array(providedGW)[bsort[ik]] for ik in xrange(NKPT)])

print "eigenDFTGW - sorted"
print eigenDFTGW

print "writing eig"
feig_out=open(seednameGW+".eig","w")
for ik in xrange(NKPT):
    for ib in xrange(NBND):
	feig_out.write(" {0:4d} {1:4d} {2:17.12f}\n".format(ib+1,ik+1,eigenDFTGW[ik,ib]))
feig_out.close()


if calcAMN:
  try:
    print "----------\n AMN   \n---------\n"
    f_amn_out=open(seednameGW+".amn","w")
    f_amn_in=open(seedname+".amn","r")
    s=f_amn_in.readline().strip()
    f_amn_out.write("{0}, sorted by quasi-particle energies (YAMBO) on {1} \n".format(s,datetime.datetime.now().isoformat()) )
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
    print "----------\n MMN  \n---------\n"

    f_mmn_out=open(os.path.join(seednameGW+".mmn"),"w")
    f_mmn_in=open(os.path.join(seedname+".mmn"),"r")

    s=f_mmn_in.readline().strip()
    f_mmn_out.write("{0}, sorted by quasi-particle energies (YAMBO) on {1} \n".format(s,datetime.datetime.now().isoformat()) )
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
    print "----------\n  {0}  \n---------".format(ext)
    
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
		    f_uXu_out.write("\n".join(x for x in (A[start:start+length].reshape(nbnd,nbnd,order='F')[BANDSORT[KPNB[ik][ib2]],:][:,BANDSORT[KPNB[ik][ib1]]]+
			np.einsum('ln,lm,l->nm',MMN[ik][ib2].conj(),MMN[ik][ib1],eigenDE[ik]) ).reshape(-1,order='F') )+"\n")
		    start+=length
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
    print "----------\n {0} OK  \n---------\n".format(ext)
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



