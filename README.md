# One-electron Diatomic Molecule Structure Calculation wf.f
The program wf.f is the combined version of 1, power.f [1] and 2, grave.f [2]. The main purpose of the program is to calculate the potential energy and the wavefunction expansion terms, which are later used to determine the couplings. 

#system requirements
The program is successfully compiled with GNU Fortran (Homebrew GCC 12.2.0) 12.2.0 on macOS Monterey version 12.6.2. It is recommended to add the -w flag when compiling, i.e.:
                                          gfortran -W wf.f 
to suppress the warnings from Fotran77. 
