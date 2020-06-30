#!/usr/bin/python
from math import *
import numpy as np
import sys
import os
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.analyzer import SpacegroupOperations
from pymatgen.core.operations import SymmOp

def Volume(a, b, c):
    v=np.array([[a[0],a[1],a[2]],
                [b[0],b[1],b[2]],
                [c[0],c[1],c[2]]])
    return (np.linalg.det(v))

def distance(x1,y1,z1,x2,y2,z2,scale,a,b,c):
    i=((x1-x2)*a[0]+(y1-y2)*b[0]+(z1-z2)*c[0])
    j=((x1-x2)*a[1]+(y1-y2)*b[1]+(z1-z2)*c[1])
    k=((x1-x2)*a[2]+(y1-y2)*b[2]+(z1-z2)*c[2])
    dist=pow((pow(i,2)+pow(j,2)+pow(k,2)),0.5)*scale
    return (dist)

def checkDistance(x,y,z,pos,R,scale,a,b,c):
    for i in pos:
        if distance(x,y,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x,y+1,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x,y-1,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x,y,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x,y,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y+1,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x,y+1,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y-1,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x,y+1,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y+1,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x,y-1,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y-1,z,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x,y-1,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y+1,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y+1,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y-1,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y+1,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y-1,z+1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y+1,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x+1,y-1,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
        if distance(x-1,y-1,z-1,float(i[0]),float(i[1]),float(i[2]),scale,a,b,c)<R:
            return (1)
    return (0)

def writeStructure(title,scale,a,b,c,at_labels,nElements,ats,ntot,pos):
    print ('%s'%(title), file=open('POSCAR','a'))
    print ('%f'%(scale), file=open('POSCAR','a'))
    print ('%f %f %f'%(a[0],a[1],a[2]),file=open('POSCAR','a'))
    print ('%f %f %f'%(b[0],b[1],b[2]),file=open('POSCAR','a'))
    print ('%f %f %f'%(c[0],c[1],c[2]),file=open('POSCAR','a'))
    if nElements==1:
        print ('%s'%(at_labels[0]), file=open('POSCAR','a'))
        print ('%d'%(ats[0]), file=open('POSCAR','a'))
    elif nElements==2:
        print ('%s %s'%(at_labels[0],at_labels[1]), file=open('POSCAR','a'))
        print ('%d %d'%(ats[0],ats[1]), file=open('POSCAR','a'))
    elif nElements==3:
        print ('%s %s %s'%(at_labels[0],at_labels[1],at_labels[2]), file=open('POSCAR','a'))
        print ('%d %d %d'%(ats[0],ats[1],ats[2]), file=open('POSCAR','a'))
    elif nElements==4:
        print ('%s %s %s %s'%(at_labels[0],at_labels[1],at_labels[2],at_labels[3]), file=open('POSCAR','a'))
        print ('%d %d %d %d'%(ats[0],ats[1],ats[2],ats[3]), file=open('POSCAR','a'))
    elif nElements==5:
        print ('%s %s %s %s %s'%(at_labels[0],at_labels[1],at_labels[2],at_labels[3],at_labels[4]), file=open('POSCAR','a'))
        print ('%d %d %d %d %d'%(ats[0],ats[1],ats[2],ats[3],ats[4]), file=open('POSCAR','a'))
    elif nElements==6:
        print ('%s %s %s %s %s %s'%(at_labels[0],at_labels[1],at_labels[2],at_labels[3],at_labels[4],at_labels[5]), file=open('POSCAR','a'))
        print ('%d %d %d %d %d %d'%(ats[0],ats[1],ats[2],ats[3],ats[4],ats[5]), file=open('POSCAR','a'))
    else:
        raise NameError('Need more space to accmmodate different kinds of elements')
    print ('Direct', file=open('POSCAR','a'))
    for m in range(ntot):
        print ('%f %f %f'%(pos[m][0],pos[m][1],pos[m][2]), file=open('POSCAR','a'))
    return()

def writeStructureLocalEnvironandAddAnAtom(x,y,z,ii,title,scale,a,b,c,at_labels,nElements,ats,ntot,pos,fft):
    lMin=12#minimun length of an axis in the expanded structure
    nx=int(ceil(lMin/distance(0,0,0,1,0,0,scale,a,b,c)))
    ny=int(ceil(lMin/distance(0,0,0,0,1,0,scale,a,b,c)))
    nz=int(ceil(lMin/distance(0,0,0,0,0,1,scale,a,b,c)))
    writeStructure(title,scale,a,b,c,at_labels,nElements,ats,ntot,pos)
    poscar = Poscar.from_file("POSCAR")
    structure = poscar.structure
    os.system('rm POSCAR*')
    fft1=np.zeros(3)
    structure.make_supercell([nx, ny, nz]);fft1[0]=fft[0]*nx;fft1[1]=fft[1]*ny;fft1[2]=fft[2]*nz#expand structure and FFT grid
    structure.to(filename="POSCAR_expanded")
    f=open('POSCAR_expanded')
    #read lattice information
    title1 = f.readline().strip()
    scale1 = f.readline()
    scale1 = float(scale1)
    a1 = f.readline().split()
    b1 = f.readline().split()
    c1 = f.readline().split()
    a1 = [float(i) for i in a1]
    b1 = [float(i) for i in b1]
    c1 = [float(i) for i in c1]
    #read atom labels and number of atoms
    at_labels1 = f.readline().split()
    nElements1=len(at_labels1)
    ats1 = f.readline().split()
    ats1 = [int(i) for i in ats1]
    ntot1 = sum(ats1)
    f.readline()
    #read positions of atoms
    pos1 = []
    for i in range(ntot1):
        d1 = f.readline().split()
        d1[0]=float(d1[0]);d1[1]=float(d1[1]);d1[2]=float(d1[2])
        pos1.append(d1)
    f.close()
    os.system('rm POSCAR*')
    rCritical=4.05#critical distance within which atoms in structures are included in the set of local environment
    nLocal=np.zeros(nElements1)#number of atoms of one specific element that is included in local environment
    nRead=0#number of atoms that have been screened
    index=[]
    for i in range(nElements1):
        for j in range(ats1[i]):
            if checkDistance(float(pos1[nRead][0]),float(pos1[nRead][1]),float(pos1[nRead][2]),[[float(x/fft1[0]),float(y/fft1[1]),float(z/fft1[2])]],rCritical,scale1,a1,b1,c1):
                index.append(nRead)
                nLocal[i]+=1
            nRead+=1
    zeroElement=[i for i, v in enumerate(nLocal) if v==0]#record elements that doesn't appear in local environment
    zeroElement.reverse()#to make sure that sequence in at_labels1 is not affected by deleting elements
    nLocalElements=nElements1-len(zeroElement)
    if len(zeroElement) != 0:
        [at_labels1.pop(i) for i in zeroElement]
        [ats1.pop(i) for i in zeroElement]
    print ('%s'%(title1), file=open('POSCAR_%s'%(str(ii)),'a'))
    print ('%f'%(scale1), file=open('POSCAR_%s'%(str(ii)),'a'))
    print ('%f %f %f'%(a1[0],a1[1],a1[2]),file=open('POSCAR_%s'%(str(ii)),'a'))
    print ('%f %f %f'%(b1[0],b1[1],b1[2]),file=open('POSCAR_%s'%(str(ii)),'a'))
    print ('%f %f %f'%(c1[0],c1[1],c1[2]),file=open('POSCAR_%s'%(str(ii)),'a'))
    if nLocalElements==1:
        print ('%s He'%(at_labels1[0]), file=open('POSCAR_%s'%(str(ii)),'a'))
        print ('%d 1'%(nLocal[0]), file=open('POSCAR_%s'%(str(ii)),'a'))
    elif nLocalElements==2:
        print ('%s %s He'%(at_labels1[0],at_labels1[1]), file=open('POSCAR_%s'%(str(ii)),'a'))
        print ('%d %d 1'%(nLocal[0],nLocal[1]), file=open('POSCAR_%s'%(str(ii)),'a'))
    elif nLocalElements==3:
        print ('%s %s %s He'%(at_labels1[0],at_labels1[1],at_labels1[2]), file=open('POSCAR_%s'%(str(ii)),'a'))
        print ('%d %d %d 1'%(nLocal[0],nLocal[1],nLocal[2]), file=open('POSCAR_%s'%(str(ii)),'a'))
    elif nLocalElements==4:
        print ('%s %s %s %s He'%(at_labels1[0],at_labels1[1],at_labels1[2],at_labels1[3]), file=open('POSCAR_%s'%(str(ii)),'a'))
        print ('%d %d %d %d 1'%(nLocal[0],nLocal[1],nLocal[2],nLocal[3]), file=open('POSCAR_%s'%(str(ii)),'a'))
    elif nLocalElements==5:
        print ('%s %s %s %s %s He'%(at_labels1[0],at_labels1[1],at_labels1[2],at_labels1[3],at_labels1[4]), file=open('POSCAR_%s'%(str(ii)),'a'))
        print ('%d %d %d %d %d 1'%(nLocal[0],nLocal[1],nLocal[2],nLocal[3],nLocal[4]), file=open('POSCAR_%s'%(str(ii)),'a'))
    elif nLocalElements==6:
        print ('%s %s %s %s %s %s He'%(at_labels[0],at_labels1[1],at_labels1[2],at_labels1[3],at_labels1[4],at_labels1[5]), file=open('POSCAR_%s'%(str(ii)),'a'))
        print ('%d %d %d %d %d %d 1'%(nLocal[0],nLocal[1],nLocal[2],nLocal[3],nLocal[4],nLocal[5]), file=open('POSCAR_%s'%(str(ii)),'a'))
    print ('Direct', file=open('POSCAR_%s'%(str(ii)),'a'))
    for m in index:
        print ('%f %f %f'%(pos1[m][0],pos1[m][1],pos1[m][2]), file=open('POSCAR_%s'%(str(ii)),'a'))
    print ('%f %f %f'%(float(x/fft1[0]),float(y/fft1[1]),float(z/fft1[2])),file=open('POSCAR_%s'%(str(ii)),'a'))
    return()

#initailzation
chg_dir='CHG'
root_dir='orgn/acetaldehyde'#'sample' for training, struct name for test
os.system('cp atom_init.json %s/'%(root_dir))
train_or_predict=False#true for train, false for predict
files=['CHG']
spacing=[0.5,0.5,0.5]#length between grid points along one axis
#files= os.listdir('./%s'%(chg_dir))#read POSCAR files of primitive cell
nStruct=0#number of structures counted
for POS in files:
    struct=POS.replace("'", '')
#    struct='aa_graphite'#name of the structure
    print (struct)
    f=open('%s/acetaldehyde/%s'%(chg_dir,struct),'r')#need for revision for training/test
    #read lattice information
    title = f.readline().strip()
    scale = f.readline()
    scale = float(scale)
    a = f.readline().split()
    b = f.readline().split()
    c = f.readline().split()
    a = [float(i) for i in a]
    b = [float(i) for i in b]
    c = [float(i) for i in c]
    V=Volume(a, b, c)#calculate volume of the unit cell
    #read atom labels and number of atoms
    at_labels = f.readline().split()
    nElements=len(at_labels)
    ats = f.readline().split()
    ats = [int(i) for i in ats]
    ntot = sum(ats)
    #ats.append(int(1))#append a hypothetical atom
    f.readline()
    #read positions of atoms
    pos = []
    for i in range(ntot):
        d = f.readline().split()
        d = [float(i) for i in d]
        pos.append(d)
    f.readline()
    #read FFT grid
    fft = f.readline().split()
    fft = [int(i) for i in fft]
    npoints=1
    for i in fft:
        npoints*=i#get total number of grids
    #read symmetry information
    writeStructure(title,scale,a,b,c,at_labels,nElements,ats,ntot,pos)
    poscar = Poscar.from_file("POSCAR")
    structure = poscar.structure
    Symmetry=SpacegroupAnalyzer(structure, symprec=0.1, angle_tolerance=5)
    symmops=Symmetry.get_symmetry_operations(cartesian=False)
    os.system('rm POSCAR*')
    #read charge density on each grid
    chg=[]
    lines=int(npoints/10)#number of lines in CHG with charge information
    if (npoints%10)!=0:
        lines+=1#to rule out possibility that the last line is not read as the last lane contain less than 10 data
    for i in range(lines):
        cd = f.readline().split()
        cd = [float(j) for j in cd]
        for k in range(len(cd)):
            chg.append(float(cd[k]/V))#charge density normalized by volume
    chg.reverse()
    f.close()

    Charge=np.zeros(shape=(fft[0],fft[1],fft[2]))
    n=0#number of grids done
    m=0#number of points included
    length_x=distance(0,0,0,1,0,0,scale,a,b,c);nx=int(round(fft[0]/int(round(length_x/spacing[0]))))#0.25 is the distance between grid points along an axis,sample one layer in nx layers in x direction
    length_y=distance(0,0,0,0,1,0,scale,a,b,c);ny=int(round(fft[1]/int(round(length_y/spacing[1]))))
    length_z=distance(0,0,0,0,0,1,scale,a,b,c);nz=int(round(fft[2]/int(round(length_z/spacing[2]))))
    if not train_or_predict:
        struct=''#faciliate csv to chg
    #transform CHGCAR to cif
    flag=False
    for k in range(0,fft[2]):
        if flag:
            print ('this material is done')
            break
        for j in range(0,fft[1]):
            if flag:
                break
            for i in range(0,fft[0]):
                if flag:
                    break
    #sample one layer in nx,ny,nz layers in x,y,z direction
                if m > 1000 and train_or_predict:
                    print ('too much data on a single structure')
                    flag=True
                    continue
                if int(i%nx)!=0 or int(j%ny)!=0 or int(k%nz)!=0:
#                print ('%d %d %d is omiited with index %d'%(i, j, k, n))
                    n+=1
                    chg.pop()
                    continue
                if Charge[i][j][k]!=0:
#                    print ('%d %d %d symmetrically equivalent'%(i, j, k))
                    n+=1
                    chg.pop()
                    ss+=1
                    if ss > 200 and train_or_predict:
                        flag=True
                        break
                    continue
                ss=0
                if checkDistance(i/fft[0],j/fft[1],k/fft[2],pos,0.01,scale,a,b,c):
#                    print ('too close to the structure')
                    n+=1
                    chg.pop()
                    continue
                charge=chg.pop()
                writeStructureLocalEnvironandAddAnAtom(i,j,k,n,title,scale,a,b,c,at_labels,nElements,ats,ntot,pos,fft)
                p = Poscar.from_file('POSCAR_%s'%(str(n)))
                w = CifWriter(p.structure, symprec=1e-6)
                w.write_file('%s/%s%08d.cif'%(root_dir,struct,n))
                os.system('rm POSCAR*') 
                print ('%s: %d %d %d finished with charge %f and index %d, %d %d %d remained' %(struct, i, j, k, charge, n,  fft[0]-i, fft[1]-j, fft[2]-k))
                print ('%s%08d,%4f' %(struct,n,charge), file=open('%s/id_prop.csv'%(root_dir),'a'))#for the convenience of recovering from .csv to CHGCAR
                n+=1
                m+=1
                Charge[i][j][k]=charge
                equ_points=[]#record equivalent points of i,j,k1
                for ops in symmops:
                    equ_points.append(ops.operate([i/fft[0],j/fft[1],k/fft[2]]))
                for equ_point in equ_points:
                    equ_a=int(round(equ_point[0]*fft[0]))
                    equ_b=int(round(equ_point[1]*fft[1]))
                    equ_c=int(round(equ_point[2]*fft[2]))
                    if equ_a >= fft[0] or equ_a<= -fft[0]:
                        equ_a -= fft[0]*int(round((equ_a/fft[0])))
#                    if equ_a < 0:
#                        equ_a += fft[0]*int(round((equ_a/fft[0])))
                    if equ_b >= fft[1] or equ_b <= -fft[1]:
                        equ_b -= fft[1]*int(round((equ_b/fft[1])))
#                    if equ_b < 0:
#                        equ_b += fft[1]*int(round((equ_b/fft[1])))
                    if equ_c >= fft[2] or equ_c <= -fft[2]:
                        equ_c -= fft[2]*int(round((equ_c/fft[2])))
#                    if equ_c < 0:
#                        equ_c += fft[2]*int(round((equ_c/fft[2])))
#                    print (equ_a,equ_b,equ_c)
                    Charge[equ_a][equ_b][equ_c]=charge
    nStruct+=1
    print ('%d structures have been counted'%(nStruct))



