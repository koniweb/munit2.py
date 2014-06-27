#!/usr/bin/env python
######################################################################
# munit.py
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
######################################################################
version="2.0"

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
    
    arg=0
    while ( (arg+1) < len(sys.argv)):
        arg=arg + 1
        if sys.argv[arg].startswith('--'):
            option = sys.argv[arg][2:]
            if   option == 'h':
                start(version)
                showhelp()
            elif option == 'in':
                inf=readopt(arg)
                arg=arg+1
            elif option == 'coo':
                file_coord=readfilename(arg)
                arg=arg+1
            elif option == 'a':
                a=[0.0,0.0,0.0]
                readvec(a,arg)
                vec[0][:]=a
                arg=arg+3
            elif option == 'b':
                b=[0.0,0.0,0.0]
                readvec(b,arg)
                vec[1][:]=b
                arg=arg+3
            elif option == 'c':
                c=[0.0,0.0,0.0]
                readvec(c,arg)
                vec[2][:]=c
                arg=arg+3
            elif option == 'off':
                off=[0.0,0.0,0.0]
                readvec(off,arg)
                vec[3][:]=off
                arg=arg+3
            elif option == 'm':
                readm(m,arg)
                arg=arg+3
            elif option == 'f':
                factor[0]=readfactor(arg)
                factor[1]=factor[0]
                factor[2]=factor[0]
                arg=arg+1
            elif option == 'fxyz':
                readvec(factor,arg)
                arg=arg+3
            elif option == 'out':
                out=readopt(arg)
                arg=arg+1
            elif option == 'q':
                quiet=1
            else:
                start(version)
                print >>sys.stderr, "option not known"
                stop()
        else:
            start(version)
            print >>sys.stderr, 'command ', sys.argv[arg], ' not known'
            stop()
    #-- START PROGRAM PROMPT ----------------------------------------------
    if quiet == 0:
        start(version)
        print >>sys.stderr  
        print >>sys.stderr,  "COORDINATES:"

    #-- READ COORDINATE FILES ---------------------------------------------
    mol=readinfo(inf,file_coord)

    # add vec and offset
    mol.set_periodicity(localvec[0],localvec[1],localvec[2],offset)

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
# read vector
def readvec(vector,argument):
    if ((argument+1+3) > len(sys.argv)):
        print >>sys.stderr, 'vectors not given completly'
        stop()
    else:
        vector[0]=float(sys.argv[argument+1])
        vector[1]=float(sys.argv[argument+2])
        vector[2]=float(sys.argv[argument+3])

# read filename
def readfilename(argument):
    if ((argument+1+1) > len(sys.argv)):
        print  >>sys.stderr, 'vector not given completly'
        stop()
    else:
        return sys.argv[argument+1]

# readm
def readm(m,argument):
    if ((argument+1+3) > len(sys.argv)):
        print  >>sys.stderr, 'm not given completly'
        stop()
    else:
        m[0]=int(sys.argv[argument+1]) #mx
        m[1]=int(sys.argv[argument+2]) #my
        m[2]=int(sys.argv[argument+3]) #mz

# read factor
def readfactor(argument):
    if ((argument+1+1) > len(sys.argv)):
        print  >>sys.stderr, 'factor not given completly'
        stop()
    else:
        return float(sys.argv[argument+1]) 
        print  >>sys.stderr, 'subrout', float(sys.argv[argument+1]) 

# reading options
def readopt(argument):
    if ((argument+1+1) > len(sys.argv)):
        print  >>sys.stderr, 'output or input type not given completly'
        stop()
    else:
        if   (sys.argv[argument+1] == 'xyz'):
            return 'xyz'
        elif (sys.argv[argument+1] == 'lammps'):
            return 'lammps'
        elif (sys.argv[argument+1] == 'pwscfin'):
            return 'pwscfin'
        elif (sys.argv[argument+1] == 'pwscfout'):
            return 'pwscfout'
        elif (sys.argv[argument+1] == 'pwscf'):
            return 'pwscf'
        else:
            print  >>sys.stderr, 'out- or input option not known'
            stop()

# printing info
def printinfo(file_coord,mol,m,factor):
    print >>sys.stderr, ("{:15s} {:<20s}".format("coordinatefile:", file_coord))
    print >>sys.stderr, ("{:15s} {:<20d}".format("natoms:",mol.natoms))
    v=mol.vec
    o=mol.offset
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
    for i in range(mol.natoms):
        print  >>sys.stderr, ("{:<6d} {:5s} {:f} {:f} {:f}".format(
                i, mol.at[i].name,
                mol.at[i].coord[0],
                mol.at[i].coord[1],
                mol.at[i].coord[2]))
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
        mol=mol.readpwscf("")       
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
    print >>sys.stderr, "#   by kweber 2012"
    print >>sys.stderr, "#   version", version
    print >>sys.stderr, "#-----------------------------"
def stop():
    print >>sys.stderr, 'While starting munit.sh an error occured. Type "munit.py --h" for help'
    sys.exit()
def showhelp():
    print >>sys.stderr, 'to invoke munit.py you have the following options:'
    print >>sys.stderr, '--h                   show this help'
    print >>sys.stderr, '--q                   quiet output'
    print >>sys.stderr, '--coo  <coord>        show coordinate file'
    print >>sys.stderr, '--a    <x> <y> <z>    x, y, z components from unitvector a'
    print >>sys.stderr, '--b    <x> <y> <z>    x, y, z components from unitvector b'
    print >>sys.stderr, '--c    <x> <y> <z>    x, y, z components from unitvector c'
    print >>sys.stderr, '--off  <x> <y> <z>    x, y, z offset of unitvectors'
    print >>sys.stderr, '--f    <factor>       stretch factor for the unitcell '
    print >>sys.stderr, '--fxyz <x> <y> <z>    stretch factor for the unitcell '
    print >>sys.stderr, '--m    <x> <y> <z>    multiplication in x, y, z direction'
    print >>sys.stderr, '--out  <option>       option for output'
    print >>sys.stderr, '                        xyz -- standard (no multiplexyz files)'
    print >>sys.stderr, '                        lammps'
    print >>sys.stderr, '                        pwscf'
    print >>sys.stderr, '--in   <option>       option for input'
    print >>sys.stderr, '                        xyz -- standard'
    print >>sys.stderr, '                        lammps'
    print >>sys.stderr, '                        pwscfin'
    print >>sys.stderr, '                        pwscfout'
    print >>sys.stderr, ''
    print >>sys.stderr, 'example:'
    print >>sys.stderr, 'munit.py --coo test.xyz --a 1 0 0 --b 0 1 0 --c 0 0 1 --m 1 1 1 '
    sys.exit()


######################################################################
# MAIN PROGRAM START
######################################################################
if __name__ == '__main__':
    coord = ''
    main()
