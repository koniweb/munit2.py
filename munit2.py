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
# 2014-07-30 Version 2.2  -- lammps input reads charge and molid     #
# 2014-08-07 Version 2.3  -- multiple xyz output                     #
# 2014-08-14              -- check for input file                    #
# 2014-11-11 Version 2.4  -- Correction of pwscfdata input reading   #
# 2014-12-23 Version 2.5  -- sort atoms and choose collection        #
# 2015-01-30 Version 2.6  -- select structures out of files          #
# 2015-02-13 Version 2.7  -- select atoms in molecules               #
######################################################################
version="2.7"

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
    inf=['xyz']
    out=['xyz',False]
    verbose=int(0);
    datapwscf=""
    structures=[]
    # vector, localvec, offset
    factor=[float(1.000),float(1.000),float(1.000)]
    localvec=[ [float(1.000),float(0.000),float(0.000)],
               [float(0.000),float(1.000),float(0.000)],
               [float(0.000),float(0.000),float(1.000)]]
    offset=    [float(0.000),float(0.000),float(0.000)]
    m=[int(0),int(0),int(0)]
    # sort variable
    sortdir=[]
    # selection variable
    selection=[]

    #-- PRINT START INFORMATION -------------------------------------------
    start(version)   

    # READING user commands
    if len(sys.argv) < 2:
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
            if len(arg)>2:  
                for i in range(2,len(arg)):
                    structures.append(int(arg[i]))
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
        # sort  datafield sortdir
        elif arg[0]=='s':
            if len(arg)==1:
                sortdir.append(0)
                sortdir.append(1)
                sortdir.append(2)
            for argi in arg[1:]: 
                if   argi=="x":sortdir.append(0)
                elif argi=="y":sortdir.append(1)
                elif argi=="z":sortdir.append(2)
                else: 
                    print >> sys.stderr, ("sorting direction \"{:s}\" not known").format(argi)
                    stop()
        # select atoms
        elif arg[0]=='sel':
            readselection(arg[1:],selection)
        # verbose
        elif arg[0]=='v':
            verbose=1
            if len(arg)>1: 
                if (arg[1]=="xyz"): verbose=2
                else: print >>sys.stderr, ("verbose option:\"{:s}\" not known").format(arg[1])
        # else
        else:
            print >>sys.stderr, ("option \"{:s}\" not known").format(arg[0])
            stop()
            
    
    #-- CHECK FOR INPUT FILE ----------------------------------------------
    # filenames
    try:
        testopen=open(file_coord,"r")
    except (IOError, OSError,NameError):
        print >>sys.stderr, "ERROR: Input file not given"
        stop()
    testopen.close()

    #-- READ COORDINATE FILES ---------------------------------------------
    mol=readinfo(inf,file_coord)

    # do selection of structures
    if len(structures)>0:
        molnew=[]
        text=""
        for struct in structures:
            if struct < len(mol): 
                molnew.append(mol[struct])
                text+=" {:d}".format(struct)
            else:
                print >>sys.stderr, "ERROR: molecule does not have {:d} structures".format(struct)
        mol=molnew
        print >> sys.stderr, "...structures {:s} were excerpted".format(text)

    # check if multiple xyz
    if ( (( not out[0]=="xyz" ) or ( out[0]=="xyz" and out[1]==False)) ):
        if (len(mol)>=1): mol=[mol[-1]]
        else: 
            print >>sys.stderr, "ERROR: Input file does not contain a structure"
            stop()        

    #-- CHECK FOR M -------------------------------------------------------
    if verbose > 0 and m == [0,0,0]:
        print >>sys.stderr,  "... multiplication automatically set to 1 1 1 \n"
        m=[1,1,1]

    #-- OUTPUT MOLECULES --------------------------------------------------
    # loop over all molecules in input file
    cntmol=0
    for moli in mol:
        cntmol+=1
        # print atomcounter
        print >>sys.stderr, ('# Molecule {:d}/{:d}').format(cntmol,len(mol))
        if verbose > 0:
            print >>sys.stderr, ('#-----------------------------------------------------------')
            
        # read data from pwscf file
        if not datapwscf=="":
            moli.setup_pwscf=moli.SETUP_PWSCF()
            moli.read_setup_pwscf(datapwscf)                                       
        
        # add vec and offset
        v=moli.vec()
        if (v[0]==[0.0,0.0,0.0] and 
            v[1]==[0.0,0.0,0.0] and 
            v[2]==[0.0,0.0,0.0]):
            moli.set_vecs(a=localvec[0],b=localvec[1],c=localvec[2])
        if moli.offset==[0.0,0.0,0.0]:
            moli.set_vecs(off=offset)        

        #-- stretch -----------------------------------------------------------
        moli.stretch(factor)
        
        #-- multiply ----------------------------------------------------------
        moli.mol_multiply(m[0],m[1],m[2])
        
        #-- print info---------------------------------------------------------
        if verbose > 0: printinfo(file_coord,moli,m,factor,sortdir,selection)
        if verbose > 1: printcoo(file_coord,moli,m,factor)
        
        # sorting
        for dir in sortdir:
            print >> sys.stderr, ("...sorting in {:d}-direction").format(dir)
            moli.sortatoms(dir,-1)

        # selection
        if len(selection>0):
            atomlist=[]
            deletelist=[]
            # build full list
            for sel in selection:
                if sel[0]==-1: sel[0]=0
                if sel[1]==-1: sel[1]=moli.natoms()-1
                for i in range(sel[0],sel[1]+1): atomlist.append(i)
            # reverse selection
            for i in range(moli.natoms()):
                if not (i in atomlist): deletelist.append(i)
            # delete atoms from deletelist
            deletelist.sort(reverse=True)
            moli.delete_atoms(deletelist)
            #print >> sys.stderr,deletelist  #DEBUG

        #-- output ------------------------------------------------------------
        output(version,out,moli)
        

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
        # XYZ
        if   (arg[1] == 'xyz'):
            if    len(arg)>2 and arg[2]=="e": return ['exyz']
            else: return ['xyz']                    
        # LAMMPS
        elif (arg[1] == 'lammps'):
            # check for additional options
            lchargein=False
            lmolin=False
            t="lammps"
            for iarg in range(2,len(arg)):
                if arg[iarg]=="c":   lchargein=True
                elif arg[iarg]=="m": lmolin=True
                elif arg[iarg]=="out": t=t+"out"
            return [t,lchargein,lmolin]
        # PWSCF
        elif (arg[1] == 'pwscf'):
            if len(arg)==3 and arg[2]=="out":
                return ['pwscfout']
            else:
                return ['pwscfin']
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
            # check for additional options
            lmultxyz=False
            for iarg in range(2,len(arg)):
                if arg[iarg]=="m":lmultxyz=True
            return ['xyz',lmultxyz]
        elif (arg[1] == 'lammps'):
            # check for additional options
            lchargeout=False
            lmolout=False
            t="lammps"
            for iarg in range(2,len(arg)):
                if arg[iarg]=="c":   lchargeout=True
                elif arg[iarg]=="m": lmolout=True
            return [t,lchargeout,lmolout]
        elif (arg[1] == 'pwscf'):
            return ['pwscf']
        else:
            print  >>sys.stderr, 'output option not known'
            stop()

# printing info
def printinfo(file_coord,mol,m,factor,sortdir,selection):
    S="  "
    #INFORMATION
    print >>sys.stderr, ("{:s}INFORMATION:").format(S)
    print >>sys.stderr, ("{:s}{:s}{:15s} {:<20s}".format(S,S,"coordinatefile:", file_coord))
    print >>sys.stderr, ("{:s}{:s}{:15s} {:<20d}".format(S,S,"natoms:",mol.natoms()))
    v=mol.vec()
    o=mol.offset()
    print >>sys.stderr, ("{:s}{:s}{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format(S,S,"vector a:", v[0][0], v[0][1], v[0][2]))
    print >>sys.stderr, ("{:s}{:s}{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format(S,S,"vector b:", v[1][0], v[1][1], v[1][2]))
    print >>sys.stderr, ("{:s}{:s}{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format(S,S,"vector c:", v[2][0], v[2][1], v[2][2]))
    print >>sys.stderr, ("{:s}{:s}{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format(S,S,"offset:", o[0], o[1], o[2]))
    print >>sys.stderr
    f=factor
    #OPTIONS
    print >>sys.stderr, ("{:s}{:15s}".format(S,"OPTIONS:"))
    print >>sys.stderr, ("{:s}{:s}{:15s} {:<15.10f} {:<15.10f} {:<15.10f}".format(S,S,"factor:", f[0], f[1], f[2]))
    print >>sys.stderr, ("{:s}{:s}{:15s} {:<5d} {:<5d} {:<5d}".format(S,S,"multiplication:", m[0],m[1],m[2]))
    # sorting
    if len(sortdir)>0:
        print >>sys.stderr, ("{:s}{:s}{:15s} {:<5d} {:<5d} {:<5d}".format(S,S,"sorting:", sortdir[0],sortdir[1],sortdir[2]))
    #else: print >>sys.stderr, ("{:s}{:s}{:15s} {:<5s}".format(S,S,"sorting:", "NONE"))
    # selection
    for sel in selection:
        print >>sys.stderr, ("{:s}{:s}{:15s} {:<5d} {:<5d}".format(S,S,"selection:", sel[0],sel[1]))
    #else: print >>sys.stderr, ("{:s}{:s}{:15s} {:<5s}".format(S,S,"selection:", "ALL"))
    print >>sys.stderr

# printing coordinates
def printcoo(file_coord,mol,m,factor):
    S="  "
    # COORDINATES
    print >>sys.stderr, ("{:s}COORDINATES:").format(S)
    for i in range(mol.natoms()):
        print  >>sys.stderr, ("{:s}{:s}{:<6d} {:5s} {:f} {:f} {:f}".format(
                S,S,
                i, mol.at()[i].type()[0],
                mol.at()[i].coord()[0],
                mol.at()[i].coord()[1],
                mol.at()[i].coord()[2]))
    print >>sys.stderr 

# READING coordinate file
def readinfo(inf,file_coord):
    filetype=inf[0]
    mol=cm.molecule()
    # XYZ
    if   (filetype=="xyz"):
        mol=mol.readxyz(file_coord)
    elif (filetype=="exyz"):
        mol=mol.readxyz(file_coord,extended=True)
    # LAMMPS
    elif (filetype=="lammps"):
        lchargein=inf[1]
        lmolin=inf[2]
        mol=mol.readlmp(file_coord,lchargein,lmolin)
    elif (filetype=="lammpsout"):
        mol=mol.readlmpcustomout(file_coord)
    # PWSCF
    elif (filetype=="pwscfin"):
        mol=mol.readpwscfin(file_coord)
    elif (filetype=="pwscfout"):
        mol=mol.readpwscfout(file_coord)       
    else:
        print >>sys.stderr, "coordinatefile type not defined"
        stop()
    return mol

# read selection
def readselection(arg,selection):
    # TODO check for erreous selections
    for argi in arg:
        argi=argi.split(":")
        # substitute empty string with -1 for beginning or end
        for i in range(len(argi)):
            if not argi[i]: argi[i]=-1
        # append to selection array
        if len(argi)==2: selection.append([int(argi[0])-1,int(argi[1])-1])

# write output file
def output(version,out,mol):
    filetype=out[0]
    if   (filetype=="xyz"):
        mol.writexyz("",extended=True)
    elif (filetype=="lammps"):
        lchargeout=out[1]
        lmolout=out[2]
        mol.writelmp("",lcharge=lchargeout,lmoltype=lmolout)
    elif (filetype=="pwscf"):
        mol.writepwscf("")       
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
    print >>sys.stderr, '--v [xyz]                verbose output'
    print >>sys.stderr, '--coo  <coord>           show coordinate file'
    print >>sys.stderr, '--a    <x>  [<y>   <z>]  x, y, z components from unitvector a'
    print >>sys.stderr, '--b   [<x>]  <y>  [<z>]  x, y, z components from unitvector b'
    print >>sys.stderr, '--c   [<x>   <y>]  <z>   x, y, z components from unitvector c'
    print >>sys.stderr, '--off  <x>   <y>   <z>   x, y, z offset of unitvectors'
    print >>sys.stderr, '--f    <a>  [<b>   <c>]  stretch factor for the unitcell '
    print >>sys.stderr, '--m    <x>   <y>   <z>   multiplication in x, y, z direction'
    print >>sys.stderr, '--s    [dir] [dir] [dir] sort atoms with directions [x,y,z]'
    print >>sys.stderr, '--sel  [a:b] [...]       select atoms a to b and [...]'
    print >>sys.stderr, '                         atomlist from 1 to natoms'
    print >>sys.stderr, '--in   <option>          option for input [xyz]'
    print >>sys.stderr, '                           xyz [e] -- optional extended xyz readin'
    print >>sys.stderr, '                           lammps [c] [m] --include charge/mid'
    print >>sys.stderr, '                           pwscf [in/out]'
    print >>sys.stderr, '--datapwscf <file>       file to read pwscf simulation setup'
    print >>sys.stderr, '--out  <option>          option for output [xyz]'
    print >>sys.stderr, '                           xyz [m] -- single or [m]ultiple xyz file'
    print >>sys.stderr, '                           lammps [c] [m] --include charge/mid'
    print >>sys.stderr, '                           pwscf'
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
