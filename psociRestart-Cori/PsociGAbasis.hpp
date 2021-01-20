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
 *   @file PsociGAbasis.hpp
 *
 */

#ifndef PSOCIGABASIS_H
#define PSOCIGABASIS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include <ga++.h>

#include "PsociTimer.hpp"
#include "PsociGAhamiltonian.hpp"

using namespace std;

const int MAX_GA_VECTOR_DIMS = 2;
/*
const double SQ_MIN_PRECISION = 0.000000001; // Used for determine if <a|b> is one or zero;
const double MIN_PRECISION = sqrt( SQ_MIN_PRECISION );
*/

const double MIN_PRECISION = 0.000001;
//const double MIN_PRECISION = 0.0000001;

//const double SQ_MIN_PRECISION = MIN_PRECISION * MIN_PRECISION;

 const double SQ_MIN_PRECISION = 0.0; // override this for now

/* Making this too small get us close to singular with concommitant issues in convergence
   Really need to move to Jacobi-Davidson or SPAM
*/
const double DAVIDSON_MIN = 0.001; // Do not want this too small

const int NUM_ORTHO_REPS = 1; // Number of repeats in MGS scheme.

const int COLS_DESIRED_UNSPECIFIED = -1;

class PsociGAbasis {
  
private:
  int64_t test;
  PsociGAhamiltonian * local_hamiltonian;  
  GA::GlobalArray * local_ga;
  GA::GlobalArray * local_Hga;  

  GA::GlobalArray * g_scratch; // Used inly in the preconditioner as a temp holder of update vector

  int vector_lo[2],vector_hi[2];
  int hvector_lo[2],hvector_hi[2];
  int vector_ndim;
  int vector_dims[2];

  int scratch_dims[2];
  int scratch_ndim;

  bool l_scatter_products; // false means to do a replicated H*v
  void fetchPseudoRandom( int number, int seed, vector<double> & buf );

//NEW BVERSION`
  int global_rows; // ( == nbasis == maxsef )
  int global_cols; // cols of input g_v: ( numSoughtRoots + maxsubspace )
  int existing_cols; // Reports number of vectors on entry - these are supposed to indicate current set of approxmations.
  int desired_cols; // The total numkber actually desired ( == numSoughtRoots) 
  int converged_cols; // The lowest set of vectors with rnorms < tolerance. These should be kept and unmodified.
  int actual_roots_sought; // The number of roots to optimize ( == numSoughtRoots )

  int present_matrixProducts;   //Number of valid Hv in local_g_Hv  starting at the bottom.
  int required_matrixProducts;  // Number of Hv to compute num = existing_cols - present_matrixProducts.
  void precondition( double shift, GA::GlobalArray * g_scratch, GA::GlobalArray * g_diag );
  


public:
  explicit PsociGAbasis( int numSoughtRoots, PsociGAhamiltonian * hamilton, GA::GlobalArray * g_array,  GA::GlobalArray * g_Harray );

//  ~PsociGAbasis();

  GA::GlobalArray * fetchHandle();   // THe set of V data
  GA::GlobalArray * fetchHvHandle(); // The Hv data 
  int CurrentExistingRoots(); 
  int CurrentNumBasis(); // maxsef
  int CurrentNumSoughtRoots();
  int CurrentMaxRoots(); // numSought + maxsubspace
  int resizeVectors( int setCols ); // Change the number of SoughtRoots adds/subtracts vectors accordingly relative to existing vecs
  int resizeVectors( int setCols, pair<int,double> & time  );
  int appendVectors( GA::GlobalArray * g_vectors );
  int appendVectors( GA::GlobalArray * g_vectors, pair<int,double> & time  );
  int getAllNorms( vector<double> & norms );
  void getAllNorms( vector<double> & norm, pair<int,double> & time );
  void normAllVectors(pair<int,double> & time );
  void normAllVectors();
  void normAppendedVectors( int start, int num );
  int checkAllOverlap( vector<pair<int,int> > & overlap );
  int checkAllOverlapGA( vector<pair<int,int> > & overlap );
  int orthonormalizeAllVectors();
  int orthonormalizeAllVectors( pair<int,double> & time);
  int orthonormalizeAppendedVectors(int start, int num );
  int orthonormalizeAppendedVectorsGA(int start, int num );
  void compactVectors();
  int computeHvProducts( int start, int num );
  int computeHvProductsSplit( int start, int num );
  int appendScratchVector();
  int CurrentReqHVProducts();
  int preconditionAndNewVectorIncore( int curroot, double shift, GA::GlobalArray * g_r, GA::GlobalArray * g_diag );
  int CurrentExistingHVProducts();
  double fetchMinPrecision();
  void setpresent_matrixProducts( int l_present_matrixProducts );
  int replaceVectors( GA::GlobalArray * g_vector, int num );
  int replaceVectorsNew( vector<vector<double> >&tvec, vector<vector<double> >&thvec, int num );

  int getAppendedNorms( vector<double> & norms, int start, int num );
  int getAppendedNormsDataParallel( vector<double> & norms, int start, int num );
  int expandAndReplaceVectors( GA::GlobalArray * g_evecs, int num );
  int expandAndReplaceVectorsNew( GA::GlobalArray * g_evecs, int num );
  void setMatrixVectorChunk( int nchunk );

  void prefetchCIMAT();
  void prefetchListCIMAT();
  void generateListCIMAT();
  void destroyLocalCIMAT();

  void generateStaggeredChunkList(int l_vectorchunk );
  void destroyStaggeredChunkList();
  void destroyGACIMAT();

  void allocateScratchVectorSpace();
  void allocateScratchVectorSpaceScatterAccumulate();
  void destroyScratchVectorSpace();
  void allocateScratchCIMATSpace();

  void setProducts( bool scatter_products );
  int resetVectorChunk(int maxmempercore );
  bool getProducts();



};

#endif //PSOCIGABASIS_H
