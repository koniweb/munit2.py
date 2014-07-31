Multiply UNITbox
================
version 2.2  

to invoke munit2.py you have the following options:  
  
option            | description
:-----------------|:-----------------------------------------
--h               | show this help
--q               | quiet output
--coo  coord      | show coordinate file
--a    x  [y   z] | x, y, z components from unitvector a
--b   [x]  y  [z] | x, y, z components from unitvector b
--c   [x   y]  z  | x, y, z components from unitvector c
--off  x   y   z  | x, y, z offset of unitvectors
--f    a  [b   c] | stretch factor for the unitcell  
--m    x   y   z  | multiplication in x, y, z direction
--in   option     |  * option for input
                  |  * xyz -- standard
                  |  * lammps [c] [m] --include charge/mid
                  |  * pwscf [in/out]
--datapwscf file  | file to read pwscf simulation setup
--out  option     | option for output
                  |  * xyz -- standard
                  |  * lammps [c] [m] --include charge/mid
                  |  * pwscf
                                                                  
example:
--------
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
munit2.py --coo test.xyz --a 1 --b 0 1 0 --c 0 0 1 --m 1 1 1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

todo list:
==========
[ ] extensive testing
