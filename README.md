# PSOCI
Parallel Spin Orbit configuration Interaction (PSOCI) program. 
PSOCI is a computationally intensive method that uses MPI and Global Arrays to manage 
largescale (parallel) data distribution and concurrancy. The enormous complexity limits 
the application of PSOCI to only the largest of computer systems. PSOCI has reportedly scaled well
50-100k nodes. It requires massive (aggregate) memory, large parallel disk space. It is somewhat related 
to the Columbus Suite of computational codes in that many iof the same shared libraries and data structures are used
and interchangable.

The PSOCI model is a relativistic, conventional CI, ab initio model, that buils orbital-based determinants 
from double-group symetry adapted basis orbitals, and performed 1+2 excitations from that space. 

The code is useful in detailed theoretical analysis of complex bonding situations such as encountered 
in Heavy Transition metal, Lanthanide or Actinide containing molecules.

This code is difficult for the non expert implement given, the need for manually setup.
