#!/usr/bin/env python
######################################################################
# munit.py                                                           #
# script to read and write coordinates in various formats            #
# to add: multiple xyz                                               #
#--------------------------------------------------------------------#
# Konstantin Weber                                                   #
# 2012-01-25 Version 1.0                                             #
# 2013-03-14 Version 1.1  -- added version                           #
# 2013-03-28 Version 1.2  -- formated output                         #
# 2013-04-23 Version 1.3  -- read exyz and multiply                  #
# 2013-11-27 Version 1.4  -- improve writing lammps files            #
# 2013-12-11 Version 1.5  -- error in pwscf input fixed              #
# 2013-12-11 Version 1.6  -- pwscf output readable                   #
# 2014-06-26 Version 2.0  -- shift to molecule                       #
# 2014-06-30 Version 2.1  -- allow to read data from pwscf           #
#                            TODO: TEST                              #
######################################################################
version="2.1"

#----------------------------------------------------------------------
# import
#----------------------------------------------------------------------
import sys
import math

#----------------------------------------------------------------------
# classes
#----------------------------------------------------------------------
import class_molecule as cm

#----------------------------------------------------------------------
# program environment
#----------------------------------------------------------------------
def main():
    #- variables declaration ------------------------------------------
    ndim=3
    out='xyz'
    inf='xyz'
    quiet=int(0);
    datapwscf=""
    # vector, localvec, offset
    factor=[float(1.000),float(1.000),float(1.000)]
    localvec=[ [float(1.000),float(0.000),float(0.000)],
               [float(0.000),float(1.000),float(0.000)],
               [float(0.000),float(0.000),float(1.000)]]
    offset=    [float(0.000),float(0.000),float(0.000)]
    m=[int(0),int(0),int(0)]

    # READING user commands
    if len(sys.argv) < 2:
        start(version)
        print >>sys.stderr, 'no user arguments'
        stop()
    
    # argument parsing
    args=[]
    A=" ".join(sys.argv).split("--")
    for i in range(1,len(A)): args.append(A[i].strip().split(" "))

    # go through list
    for arg in args:
        if   arg[0]=="h":
            start(version)
            showhelp()
        # input files
        elif arg[0]=='coo':
            file_coord=readfilename(arg)
        elif arg[0]=='datapwscf':
            datapwscf=readfilename(arg)
        # file types
        elif arg[0]=='in':
            inf=readin(arg)
        elif arg[0]=='out':
            out=readout(arg)
        # vectors
        elif arg[0]=='a':
            localvec[0][:]=readvec(arg,localvec[0][:],'a',0)
        elif arg[0]=='b':
            localvec[1][:]=readvec(arg,localvec[1][:],'b',1)
        elif arg[0]=='c':
            localvec[2][:]=readvec(arg,localvec[2][:],'c',2)
        elif arg[0]=='off':
            offset=readvec(arg,offset,'offset')
        # stretch and multiply options
        elif arg[0]=='f':
            factor=readfactor(arg)
        elif arg[0]=='m':
            m=readvecint(arg,m,'m')
        # quiet
        elif arg[0]=='q':
            quiet=1
        # else
        else:
            start(version)
            print >>sys.stderr, "option ",arg[0]," not known"
            stop()

    #-- START PROGRAM PROMPT ----------------------------------------------
    if quiet == 0:
        start(version)
        print >>sys.stderr  
        print >>sys.stderr,  "COORDINATES:"

    #-- READ COORDINATE FILES ---------------------------------------------
    mol=readinfo(inf,file_coord)
    
    # read data from pwscf file
    if not datapwscf=="":
        mol.setup_pwscf=mol.SETUP_PWSCF()
        mol.read_setup_pwscf(datapwscf)                                       

    # add vec and offset
    v=mol.vec()
    if (v[0]==[0.0,0.0,0.0] and 
        v[1]==[0.0,0.0,0.0] and 
        v[2]==[0.0,0.0,0.0]):
        mol.set_vecs(a=localvec[0],b=localvec[1],c=localvec[2])
    if mol.offset==[0.0,0.0,0.0]:
        mol.set_vecs(off=offset)

    #-- print info---------------------------------------------------------
    if quiet == 0: printinfo(file_coord,mol,m,factor)

    #-- stretch -----------------------------------------------------------
    mol.stretch(factor)
    
    #-- multiply ----------------------------------------------------------
    mol.mol_multiply(m[0],m[1],m[2])

    #-- output ------------------------------------------------------------
    output(version,out,mol)
    

#----------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------
# read vector either all directions or only in direction dir
def readvec(arg,vector,name="",dir=-1):
    if len(arg)==2: 
        vector=[0.0,0.0,0.0]
        vector[dir]=float(arg[1])
    elif len(arg)==4: vector=[float(arg[1]),float(arg[2]),float(arg[3])]
    else:
        print >>sys.stderr, 'vector ',name,' not given completly'
        stop()
    return vector
def readvecint(arg,vector,name=""):
    if len(arg) < 4:
        print >>sys.stderr, 'vector ',name,' not given completly'
        stop()
    else:
        vector[0]=int(arg[1])
        vector[1]=int(arg[2])
        vector[2]=int(arg[3])
    return vector

# read filename
def readfilename(argument):
    if len(argument) < 2:
        print  >>sys.stderr, 'file not given'
        stop()
    return argument[1]

# read factor
def readfactor(arg):
    if   len(arg) == 2:
        factor=[float(arg[1]),float(arg[1]),float(arg[1])]
    elif len(arg) == 4:
        factor=[float(arg[1]),float(arg[2]),float(arg[3])]
    else:
        print  >>sys.stderr, 'factor not given completly'
        stop()
    return factor

# reading options
def readin(arg):
    if  len(arg) < 2:
        print  >>sys.stderr, 'input type not given completly'
        stop()
    else:
        if   (arg[1] == 'xyz'):
            return 'xyz'
        elif (arg[1] == 'lammps'):
            return 'lammps'
        elif (arg[1] == 'pwscfin'):
            return 'pwscfin'
        elif (arg[1] == 'pwscfout'):
            return 'pwscfout'
        else:
            print  >>sys.stderr, 'input option not known'
            stop()
        

# reading options
def readout(arg):
    if  len(arg)< 2:
        print  >>sys.stderr, 'output type not given completly'
        stop()
    else:
        if   (arg[1] == 'xyz'):
            return 'xyz'
        elif (arg[1] == 'lammps'):
            return 'lammps'
        elif (arg[1] == 'pwscf'):
            return 'pwscf'
        else:
            print  >>sys.stderr, 'output option not known'
            stop()

# printing info
def printinfo(file_coord,mol,m,factor):
    print >>sys.stderr, ("{:15s} {:<20s}".format("coordinatefile:", file_coord))
    print >>sys.stderr, ("{:15s} {:<20d}".format("natoms:",mol.natoms()))
    v=mol.vec()
    o=mol.offset()
    print >>sys.stderr, ("{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format("vector a:", v[0][0], v[0][1], v[0][2]))
    print >>sys.stderr, ("{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format("vector b:", v[1][0], v[1][1], v[1][2]))
    print >>sys.stderr, ("{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format("vector c:", v[2][0], v[2][1], v[2][2]))
    print >>sys.stderr, ("{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format("offset:", o[0], o[1], o[2]))
    print >>sys.stderr
    f=factor
    print >>sys.stderr, ("{:15s}".format("OPTIONS:"))
    print >>sys.stderr, ("{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format("factor:", f[0], f[1], f[2]))
    print >>sys.stderr, ("{:15s} {:<5d} {:<5d} {:<5d}".format("multiplication:", m[0],m[1],m[2]))
    print >>sys.stderr
    # print atoms
    print >>sys.stderr, "coordinates:"
    for i in range(mol.natoms()):
        print  >>sys.stderr, ("{:<6d} {:5s} {:f} {:f} {:f}".format(
                i, mol.at()[i].type()[0],
                mol.at()[i].coord()[0],
                mol.at()[i].coord()[1],
                mol.at()[i].coord()[2]))
    print >>sys.stderr 

# READING coordinate file
def readinfo(inf,file_coord):
    mol=cm.molecule()
    if   (inf=="xyz"):
        mol=mol.readxyz(file_coord)
    elif (inf=="lammps"):
        mol=mol.readlmp(file_coord)
    elif (inf=="pwscfin"):
        mol=mol.readpwscfin(file_coord)
    elif (inf=="pwscfout"):
        mol=mol.readpwscfout(file_coord)       
    else:
        print >>sys.stderr, "coordinatefile type not defined"
        stop()
    return mol[-1]

# write output file
def output(version,out,mol):
    if   (out=="xyz"):
        mol=mol.writexyz("")
    elif (out=="lammps"):
        mol=mol.writelmp("")
    elif (out=="pwscf"):
        mol=mol.writepwscf("")       
    else:
        print >>sys.stderr, "output file type not defined"
        stop()
    return 

#----------------------------------------------------------------------
# help
#----------------------------------------------------------------------
def start(version):
    print >>sys.stderr, "#-----------------------------"
    print >>sys.stderr, "#   Multiply UNITbox"
    print >>sys.stderr, "#   by kweber 2014"
    print >>sys.stderr, "#   version", version
    print >>sys.stderr, "#-----------------------------"
def stop():
    print >>sys.stderr, 'While starting munit2.sh an error occured. Type "munit2.py --h" for help'
    sys.exit()
def showhelp():
    print >>sys.stderr, 'to invoke munit2.py you have the following options:'
    print >>sys.stderr, '--h                      show this help'
    print >>sys.stderr, '--q                      quiet output'
    print >>sys.stderr, '--coo  <coord>           show coordinate file'
    print >>sys.stderr, '--a    <x>  [<y>   <z>]  x, y, z components from unitvector a'
    print >>sys.stderr, '--b   [<x>]  <y>  [<z>]  x, y, z components from unitvector b'
    print >>sys.stderr, '--c   [<x>   <y>]  <z>   x, y, z components from unitvector c'
    print >>sys.stderr, '--off  <x>   <y>   <z>   x, y, z offset of unitvectors'
    print >>sys.stderr, '--f    <a>  [<b>   <c>]  stretch factor for the unitcell '
    print >>sys.stderr, '--m    <x>   <y>   <z>   multiplication in x, y, z direction'
    print >>sys.stderr, '--out  <option>          option for output'
    print >>sys.stderr, '                           xyz -- standard'
    print >>sys.stderr, '                           lammps'
    print >>sys.stderr, '                           pwscf'
    print >>sys.stderr, '--in   <option>          option for input'
    print >>sys.stderr, '                           xyz -- standard'
    print >>sys.stderr, '                           lammps'
    print >>sys.stderr, '                           pwscfin'
    print >>sys.stderr, '                           pwscfout'
    print >>sys.stderr, '--datapwscf <file>       file to read pwscf simulation setup'
    print >>sys.stderr, ''
    print >>sys.stderr, 'example:'
    print >>sys.stderr, 'munit2.py --coo test.xyz --a 1 --b 0 1 0 --c 0 0 1 --m 1 1 1 '
    sys.exit()


######################################################################
# MAIN PROGRAM START
######################################################################
if __name__ == '__main__':
    coord = ''
    main()
