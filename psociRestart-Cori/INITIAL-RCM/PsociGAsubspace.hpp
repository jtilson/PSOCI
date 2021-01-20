/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociGAbasis:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociGAsubspace.hpp
 *
 */

#ifndef PSOCIGASUBSPACE_H
#define PSOCIGASUBSPACE_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include <ga++.h>

#include "PsociGAbasis.hpp"
#include "PsociGAhamiltonian.hpp"

using namespace std;


class PsociGAsubspace {
  
private:
   PsociGAbasis * local_vector_set;
   int local_maxroots, local_numRoots, local_nsub, local_nbasis, local_maxsubspace;
   GA::GlobalArray * local_g_v;
   GA::GlobalArray * local_g_Hv;
   GA::GlobalArray * local_g_diag;



// Set of pointers for subspace (distributed) arrays

   GA::GlobalArray *g_s;
   GA::GlobalArray *g_vHv;
   GA::GlobalArray *g_evecs;
   GA::GlobalArray *g_evecs_scratch;
   GA::GlobalArray *g_diag;
   GA::GlobalArray *g_evals; 
   GA::GlobalArray *g_temp;
   GA::GlobalArray *g_temp2;

   GA::GlobalArray *g_rnorms; // Array of rnorms

   int existing_subspace_size;
   inline void daxpyLocal( int width, double value, double * tvec, double * hvec);
   inline int localPatchFetch( int * lo, int * hi, int maxvalue );


   int local_vHvlo[2], local_vHvhi[2]; // Used to save time in diagSubspaceGOP 
   int vector_lo[2], vector_hi[2];
   int hvector_lo[2], hvector_hi[2];


public:
//   PsociGAsubspace( PsociGAbasis * vector_set, GA::GlobalArray * g_Hv, GA::GlobalArray * g_diag  );
   PsociGAsubspace( PsociGAbasis * vector_set, GA::GlobalArray * g_diag, GA::GlobalArray * g_rnorms);
//   ~PsociGAsubspace();
   void setMaxSubspace( int maxsub );
   void allocateSubspace();
   void destroySubspace();
   int numRoots();
   int numVectorBasis();
   int numCIBasis();
   int maxSubspace();
   int fetchMaxSubspace();
   int fetchExistingSubspace();
   void generateSubspace();
   void generateIncrementalSubspace(int num);
   void diagSubspace();
   void diagSubspaceGA();
   void diagSubspaceGOP();
   int computeRnorms( vector<double> & evals, vector<double> & rnorms );
   int computeRnormsGA( vector<double> & evals,vector<double> & rnorms );
   int brdcstEvals( vector<double> & evals );
   int fetchEigenvectorRnorms( vector<double> & rnorms );
   void printVectorContents( vector<double> & rnorms );
   int expandEvecs( int num );
   int expandEvecsGA();
   GA::GlobalArray * fetchEvecsHandle();
   int fetchExistingSubspaceSize();
   int computeRnormSingle( int root, vector<double> & evals, double & rnormvec );




};

#endif //PSOCIGASUBSPACE_H
