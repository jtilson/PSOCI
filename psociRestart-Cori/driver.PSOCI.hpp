/**************************************************************************************
 * Copyright (c) 2011 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license

 * New implementation of PSOCI:

 Classes: 

 Description: 
 
 History:

**************************************************************************************/
/**
 *   @file driver.GAhamiltonian.Restart.C
 *
 */

#ifndef DRIVER_PSOCI
#define DRIVER_PSOCI

#include "PsociGAbasis.hpp"

  const string PSOCI_VERSION ="2.0";

// How many eigenvectors to construct 

   int numRoots = 40;

// How many roots should NOs be computer for?
// If set to <=0 then set = numRoots;

   int numNOroots = 1;

// Subspace optimization parameters

   int subspace_min_mult = 2;
   int subspace_max_mult = 4;
   int subspace_iterations = 1000;

   double coarse_rtol = 0.001;
//double coarse_rtol = 0.0001;
   double final_rtol = 1.0 * MIN_PRECISION; // From PsociGAbasis: Tighter than this causes linear dependencies

   bool showevals = false; // More of a debugging flag
 
// Basis filename information

   string moints_filename = "moints";
   string vector_filename = "vector_test.dat";
   string configs_filename = "configs.test";

// additional names not intended for stdin processing

   string DRA_CIfilename  = "DRA_cimat";
   string DRA_ICIfilename = "DRA_icicol";
   string DRA_NUMfilename = "DRA_number";
   string DRA_DIAGfilename = "DRA_diag";

// For supplemental only
   string DRA_CIfilename_supp  = "DRA_cimat_supp";
   string DRA_ICIfilename_supp = "DRA_icicol_supp";
   string DRA_NUMfilename_supp = "DRA_number_supp";

// Performance enhancements to matrix-vector products
// The effort grows linearly with value. But replicated memory decreases linearly
// Only used for -DSQUAREHAMILTONIAN

   int l_vectorchunk = 1;

// Option to adjust the NXTVAL chunking. Do not go to large

   int nxtval_chunk = 1;

// restart information - not fully implemented yet
// will require a valid vector_filename when implemented
 
   bool vector_restart = false;
   bool hamiltonian_restart = false;
   bool hamiltonian_store = false;

// An ugly parameter - the maximum sparsity anticipate for the hamiltonian
 
   int maxsparse = 3000;

   int supp_maxsparse = 8000; //Only used for driver.PSOCI.splittest.x
   int percent_maxwidth = 10; // Percentage of the maxsef used for supp storage  

   int l_granularity = 1; // Default width for 2D updates
   int l_numJchunks = 1; // Cache the entire J orbmap data object for H/CIdensity build
   
// Choose replicated bufV method for metrixProducts
   bool scatter_products = false;
//////////////////////////////////////////////////////////////
// Runtime details that should rarely be addressed by the user
// But access makes performance tesating, etc easier


/* Set up Global Array fabric
   heap and stack are per-core quantities - for collective operations
   LOCAL MEMORY - per core -Allocate in terms of DOUBLES for now
   Generally at most 25% memory/core. We want to leave 50% for O/S activities and 25% for C++ dynamic container allocation
*/

// Local per node data
  const unsigned long heap = 128L * 1024L * 1024L, stack = 128L * 1024L * 1024L ;

// GA object constraints
  const unsigned long maxGAMemPerCore = 2.0L * 1024L * 1024L * 1024L;

#endif DRIVER_PSOCI
