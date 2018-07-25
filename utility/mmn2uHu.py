#!/usr/bin/env python2
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
# Python2 script to find the indeces of a coarse mesh into a finer mesh
# provided they are commensurate.
#
# Written by Stepan Tsirkin (University of the Basque Country)
#

import numpy as np
import os
import datetime,sys
from scipy.io import FortranFile


if len(sys.argv)<2:
    print '''This utility calculates the matrices .uHu and/or .uIu from the .mmn matrices, and also reduces the number of bands in .amn, .mmn, .eig  and .spn files
	 Example of usage:
	  
		mmn2uHu.py seedname NBout=10 NBsum=100,200  targets=mmn,uHu formatted=uHu
		
	 only the first parameter is mandatory.
		NBout -- the number of bands in the output files. Default : all bands
		NBsum -- the number of bands in the summation. (one may specify several numbers, usefull to test convergence with the number of bands). Default:all bands
		input  -- path to the input files. Default: ./  
		output -- path to the output files 
		targets : files to write (amn,mmn,spn,uHu,uIu,eig) (default: amn,mmn,eig,uHu)
		formatted : files to write as formatted: uHu,uIu,spn,spn_in,spn_out,all (default: none)
		'''
    exit()
    
    
writeAMN=True
writeMMN=True
writeUHU=True
writeEIG=True
writeUIU=False
writeSPN=False

uHu_formatted     = False
uIu_formatted     = False
spn_formatted_out = False
spn_formatted_in  = False

PREFIX=sys.argv[1]

NB_out_list=[None]
NB_sum_list=[None]
INPUTDIR="./"
OUTDIR="reduced"




for arg in sys.argv[2:]:
    arg=arg.split("=")
    if arg[0]=="NBout": NB_out_list=[int(s) for s in arg[1].split(',')]
    if arg[0]=="NBsum": NB_sum_list=[int(s) for s in arg[1].split(',')]
    if arg[0]=="input": INPUTDIR=arg[1]
    if arg[0]=="output": OUTDIR=arg[1]
    if arg[0]=="targets":
	tarlist=arg[1].split(",")
	writeAMN="amn" in tarlist
	writeEIG="eig" in tarlist
	writeMMN="mmn" in tarlist
        writeUHU="uHu" in tarlist
	writeUIU="uIu" in tarlist
        writeSPN="spn" in tarlist
    if arg[0]=="formatted":
	tarlist=arg[1].split(",")
        uHu_formatted     =any(x in tarlist for x in ("uHu","all"))
	uIu_formatted     =any(x in tarlist for x in ("uIu","all"))
        spn_formatted_out=any(x in tarlist for x in ("spn","spn_out","all"))
        spn_formatted_in =any(x in tarlist for x in ("spn","spn_in","all"))




AMNrd=False
MMNrd=False
EIGrd=False

if not (writeEIG or writeUHU): EIGrd=True


print "----------\n MMN  read\n---------\n"

if not MMNrd:
	f_mmn_in=open(os.path.join(INPUTDIR,PREFIX+".mmn"),"r")
	MMNhead=f_mmn_in.readline().strip()
	s=f_mmn_in.readline()
	NB_in,NK,NNB=np.array(s.split(),dtype=int)
	MMN=[]
	MMNheadstrings=[]
	if writeMMN or writeUHU or writeUIU:
	  for ik in xrange(NK):
	    print ik
	    MMN.append([])
	    MMNheadstrings.append([])
	    for ib in xrange(NNB):
		s=f_mmn_in.readline()
		MMNheadstrings[ik].append(s)
		ik1,ik2=(int(i)-1 for i in s.split()[:2])
		tmp=np.array([[f_mmn_in.readline().split() for n in xrange(NB_in)] for m in xrange(NB_in) ],dtype=float)
		tmp=(tmp[:,:,0]+1j*tmp[:,:,1])
		MMN[ik].append(tmp)
	MMNrd=True

print "----------\n MMN  read - OK\n---------\n"

for NB_out in NB_out_list:
    if NB_out==None: NB_out=NB_in
    RESDIR="{0}_NB={1}".format(OUTDIR,NB_out)
    try:
	os.mkdir(RESDIR)
    except Exception as ex:
	print ex


    if writeMMN:
	f_mmn_out=open(os.path.join(RESDIR,PREFIX+".mmn"),"w")
	f_mmn_out.write("{0}, reduced to {2} bands {1} \n".format(MMNhead,datetime.datetime.now().isoformat(),NB_out) )
	f_mmn_out.write("  {0:10d}  {1:10d}  {2:10d}\n".format(NB_out,NK,NNB) )
        for ik in xrange(NK):
	    print ik
	    for ib in xrange(NNB):
		f_mmn_out.write(MMNheadstrings[ik][ib])
		for m in xrange(NB_out): 
	    	    for n in xrange(NB_out):
			f_mmn_out.write( "  {0:16.12f}  {1:16.12f}\n".format(MMN[ik][ib][m,n].real,MMN[ik][ib][m,n].imag) )
	f_mmn_out.close()
    print "----------\n MMN OK  \n---------\n"

    if not EIGrd:
	EIG=np.loadtxt(os.path.join(INPUTDIR,PREFIX+".eig"),usecols=(2,)).reshape((NK,NB_in),order='C')
	EIGrd=True
    
    if writeEIG:
	feig_out=open(os.path.join(RESDIR,PREFIX+".eig"),"w")
	for ik in xrange(NK):
	    for ib in xrange(NB_out):
    		feig_out.write(" {0:4d} {1:4d} {2:17.12f}\n".format(ib+1,ik+1,EIG[ik,ib]))
	feig_out.close()


    print "----------\n AMN   \n---------\n"

    
    if writeAMN:
      if not AMNrd:
	f_amn_in=open(os.path.join(INPUTDIR,PREFIX+".amn"),"r")
	head_AMN=f_amn_in.readline().strip()
	s=f_amn_in.readline()
	nb,nk,npr=np.array(s.split(),dtype=int)
	assert nb==NB_in
	assert nk==NK
	AMN=np.loadtxt(f_amn_in,dtype=float)[:,3:5]
	print "AMN size=",AMN.shape
	print nb,npr,nk
	AMN=np.reshape(AMN[:,0]+AMN[:,1]*1j,(nb,npr,nk) ,order='F')
	AMNrd=True
	f_amn_in.close()
	
    
	f_amn_out=open(os.path.join(RESDIR,PREFIX+".amn"),"w")
	f_amn_out.write("{0}, reduced to {2} bands {1} \n".format(head_AMN,datetime.datetime.now().isoformat(),NB_out) )
	f_amn_out.write("  {0:10d}  {1:10d}  {2:10d}\n".format(NB_out,NK,npr) )
	for ik in xrange(nk):
	    amn=AMN[:NB_out,:,ik]
	    for ipr in xrange(npr):
		f_amn_out.write(
		    "".join(" {0:4d} {1:4d} {2:4d}  {3:16.12f}  {4:16.12f}\n".format(ib+1,ipr+1,ik+1,amn[ib,ipr].real,amn[ib,ipr].imag)
			for ib in xrange(NB_out)))
	f_amn_out.close()
    print "----------\n AMN  - OK \n---------\n"


    UXUlist=[]
    if writeUHU:UXUlist.append(("uHu",uHu_formatted))
    if writeUIU:UXUlist.append(("uIu",uIu_formatted))
    print UXUlist
    if len(UXUlist)>0:
      for NB_sum in NB_sum_list:
	if NB_sum==None: NB_sum=NB_in
	for UXU in UXUlist:

	    print "----------\n  {1}  NBsum={0} \n---------".format(NB_sum,UXU[0])
	    formatted =UXU[1]

	    header=[" "]*60
	    head="{3} from mmn red to {1} sum {2} bnd {0} ".format(datetime.datetime.now().isoformat(),NB_out,NB_sum,UXU[0]) 
	    header[:len(head)]=head
	    header=header[:60]
	    print header
	    print head
	    print len(header)
	    print "".join(header)
	    if formatted:
		f_uXu_out = open(os.path.join(RESDIR,PREFIX+"_nbs={0:d}.{1}".format(NB_sum,UXU[0])), 'w')
		f_uXu_out.write("".join(header)+"\n")
		f_uXu_out.write("{0}   {1}   {2} \n".format(NB_out,NK,NNB))
	    else:
		f_uXu_out = FortranFile(os.path.join(RESDIR,PREFIX+"_nbs={0:d}.{1}".format(NB_sum,UXU[0])), 'w')
		f_uXu_out.write_record(header)
		f_uXu_out.write_record(np.array([NB_out,NK,NNB],dtype=np.int32))

	    for ik in xrange(NK):
		for ib2 in xrange(NNB):
		    for ib1 in xrange(NNB):
			eig_dum=None
			if UXU[0]=="uHu":
			    eig_dum=EIG[ik]
			if UXU[0]=="uIu":
			    eig_dum=np.ones(NB_sum)
			A=np.einsum('ml,nl,l->mn',MMN[ik][ib1][:NB_out,:NB_sum].conj(),MMN[ik][ib2][:NB_out,:NB_sum],eig_dum[:NB_sum]).reshape(-1,order='C')
			if(formatted):
			    f_uXu_out.write("".join("{0:20.10e}   {1:20.10e}\n".format(a.real,a.imag)  for a in A) )
			else:
			    f_uXu_out.write_record(A)
	    print "----------\n {0} OK  \n---------\n".format(UXU[0])
	    f_uXu_out.close()


    if writeSPN:
	
	print "----------\n SPN  \n---------\n"

	if spn_formatted_in:
	    f_spn_in = open(os.path.join(INPUTDIR,PREFIX+".spn"), 'r')
	    SPNheader=f_spn_in.readline().strip()
	    nbnd,NK=(int(x) for x in f_spn_in.readline().split())
	else:
	    f_spn_in = FortranFile(os.path.join(INPUTDIR,PREFIX+".spn"), 'r')
	    SPNheader=f_spn_in.read_record(dtype='c')
	    nbnd,NK=f_spn_in.read_record(dtype=np.int32)
	    SPNheader="".join(SPNheader)


	print SPNheader
	    
	assert nbnd==NB_in

	indm,indn=np.tril_indices(NB_in)
	indmQP,indnQP=np.tril_indices(NB_out)

	if spn_formatted_out:
	    f_spn_out = open(os.path.join(RESDIR,PREFIX+".spn"), 'w')
	    f_spn_out.write(SPNheader+"\n")
	    f_spn_out.write("{0}  {1}\n".format(NB_out,NK))
	else:
	    f_spn_out = FortranFile(os.path.join(RESDIR,PREFIX+".spn"), 'w')
	    f_spn_out.write_record(SPNheader)
	    f_spn_out.write_record(np.array([NB_out,NK],dtype=np.int32))


	for ik in xrange(NK):
	    A=np.zeros((3,nbnd,nbnd),dtype=np.complex)
	    if spn_formatted_in:
		tmp=np.array( [f_spn_in.readline().split() for i in xrange (3*nbnd*(nbnd+1)/2)  ],dtype=float)
		tmp=tmp[:,0]+1.j*tmp[:,1]
	    else:
		tmp=f_spn_in.read_record(dtype=np.complex)
	    A[:,indn,indm]=tmp.reshape(3,nbnd*(nbnd+1)/2,order='F')
	    check=np.einsum('ijj->',np.abs(A.imag))
	    A[:,indm,indn]=A[:,indn,indm].conj()
	    A=A[:,:NB_in,:NB_in][:,indnQP,indmQP].reshape(-1,order='F')
	    if check> 1e-10:
		raise RuntimeError ( "REAL DIAG CHECK FAILED : {0}".format(check) )
	    if spn_formatted_out:
		f_spn_out.write("".join("{0:26.16e}  {1:26.16e}\n".format(x.real,x.imag) for x in A) )
	    else:
		f_spn_out.write_record(A)
	print "----------\n SPN OK  \n---------\n"


