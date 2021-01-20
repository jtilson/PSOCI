# PSOCI
Parallel Spin Orbit Configuration Interaction (PSOCI) program. 
PSOCI is a computationally intensive method that uses MPI and Global Arrays to manage 
largescale (parallel) data distribution and concurrancy. The enormous complexity limits 
the application of PSOCI to only the largest of computer systems. PSOCI has reportedly scaled well
50-100k nodes. It requires massive (aggregate) memory and large parallel disk space. It is somewhat related 
to the Columbus Suite of computational codes in that many of the same shared libraries and data structures are used
and interchangable.

The PSOCI model is a relativistic, conventional CI, ab initio model, that builds orbital-based determinants 
from double-group symmetry adapted basis orbitals. In typical usage, a Full-CI in the valence space is specified followed by 1+2 excitations.
The single excitations provide reasonable energy values for spin-orbit fine splittings.

The code is useful in detailed theoretical analysis of complex bonding situations such as encountered 
in Heavy Transition metal, Lanthanide or Actinide containing molecules.

This code is difficult for the non expert implement given, the need for manually setup.
