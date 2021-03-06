Multiply UNITbox
================
version 2.3  

to invoke munit2.py you have the following options:  
  
option            | description
:-----------------|:-----------------------------------------
--h               | show this help
--v               | verbose output
--coo  coord      | show coordinate file
--a    x  [y   z] | x, y, z components from unitvector a
--b   [x]  y  [z] | x, y, z components from unitvector b
--c   [x   y]  z  | x, y, z components from unitvector c
--off  x   y   z  | x, y, z offset of unitvectors
--f    a  [b   c] | stretch factor for the unitcell  
--m    x   y   z  | multiplication in x, y, z direction
--s    [dir] [...] |sort atoms with directions [x,y,z]
--sel  [a:b] [...] |select atoms a to b and [...]
--in   option     |  * option for input
standard|  * xyz [e] --include extended xyz data
                  |  * lammps [c] [m] --include charge/mid
                  |  * pwscf [in/out]
--datapwscf file  | file to read pwscf simulation setup
--out  option     | option for output
standard                  |  * xyz [e] -- standard
                  |  * lammps [c] [m] --include charge/mid
                  |  * pwscf
                                                                  
example:
--------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
munit2.py --coo test.xyz --a 1 --b 0 1 0 --c 0 0 1 --m 1 1 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

todo list:
==========
- [X] read in extended xyz especially unitvectors  
- [X] add read in of pwscf out for multiple structures  
- [X] write extended xyz with unitvectors  
- [X] read in extended xyz files  
