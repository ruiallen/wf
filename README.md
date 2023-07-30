# One-electron Diatomic Molecule Structure Calculation wf.f
The program wf.f is the combined version of 1, power.f [1] and 2, grave.f [2]. The main purpose of the program is to calculate the potential energy and the wavefunction expansion terms, which are later used to determine the couplings. 

# system requirements
The program is successfully compiled with GNU Fortran (Homebrew GCC 12.2.0) 12.2.0 on macOS Monterey version 12.6.2. It is recommended to add the -w flag when compiling, i.e.:
                                          gfortran -W wf.f 
to suppress the warnings from Fotran77. 

Please make sure there is a directory WF to store the outputs. I create this folder manually, and the code will not check/create one for you. 

# Input Parameters
The input parameters are stored in file WF.DAT:
First Line
  QU: Nuclear Charge Ratio. For instance, QU = 1 for HH+, QU = 2 for HeH+, etc.
  N,L,M: Molecular Quantum Numbers.
All other parameters are irrelevant to the current work, and leave them as is.

Second Line
  Rstart: must be 0.
  IAC: Not sure of their functions. I change them to 15, 0, 100 and they don't seem to have any effects. 
  
Third Line
  NPT: The number of points, should be less than 300, otherwise the later subroutines would fail. I am working on increasing its capacity.
  RINCR: Step size in R. 


# Outpus:
All outputs are stored in the directory .../WF/QU/ where QU is an integer number like 1. 
Consider the state 3PP for HH+ system, then N = 3, L = 1 and M = 1. QU = 1.
Five files will be produced in /WF/1/:
  fort.311: Expansion terms that are used in later calculations.
  full_311.dat: Contains all information: potential, energy, and separation constants. ICVGY ICVGX are not important.
  pot_311.dat: Potential Energy Only.
  rvalue.dat: All R grids.
  IOEDM: Plain text version of full_311.dat. Not important, used in intermediate steps. 
  

# References
[1] Power, J. D., “Fixed nuclei two-centre problem in quantum mechanics”, Philosophical Transactions of the Royal Society of London. Series A, Mathematical and Physical Sciences 274, 663–697 (1973).
[2] Salin, A., “Calculation of wave-functions and collision matrix elements for one-electron diatomic molecules”, Computer Physics Communications 14, 121–132 (1978).



