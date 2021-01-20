//TODO incorporate proven weaker bounds checking. Also significant performance analysis must be performed

/**************************************************************************************

* Copyright (c) 2010, 2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: This method takes the current set of input vectors ( vector_set ), the current
              set of H*v ( tbd ) and the current preconditioner and computes the subspace eigenvectors.
              New eigenvectors are computed based on the rnorms. 

              PsociGAbasis handles vector inits, restarts, and IO. Generally not used here, however.
              PsociGAhamitonian handles DRA reads and computations if needed.

              I did not want this code to actually perform iterations and convergence estimates.
              it simply solves the subspace problem and computes all rnorms. Then it generates a 
              guess for each vector where rnornm > tolerance. The number of these kept is a user choice.

              vector set compression is NOT performed here. 
 
              The subspace always computes the subspace in Existingroots dimensions

              These are generally Collective calls. 
 History:

**************************************************************************************/
/**
 *   @file PsociGAsubspace.C
 *
 */

#include "PsociTimer.hpp"
#include "PsociGAsubspace.hpp"
#include "PsociBlasLapack.hpp"

#include<iomanip>

using namespace std;

/* 
   provide g_Hv since we know all the vector going in anyway this is just subspace processing
   provide precond diag terms
*/

PsociGAsubspace::PsociGAsubspace( PsociGAbasis * vector_set, GA::GlobalArray * g_diag, GA::GlobalArray * in_rnorms )
{
  local_vector_set = vector_set;
  local_g_v = vector_set->fetchHandle();      // set of vectors Sought+current subspace
  local_g_Hv =  vector_set->fetchHvHandle();  // set of vector products Sought+current subspace
  
  local_g_diag = g_diag;
  
  g_rnorms = in_rnorms;

  existing_subspace_size = 0; // this keeps track of the existing vHv.vSv entries
  
  char handleVecMessage[] = "PsociGAsubspace: Invalid Vectors";
  local_g_v->checkHandle( handleVecMessage );
  
  char handleCIMessage[] = "PsociGAsubspace: Invalid h_Hv mat";
  local_g_Hv->checkHandle( handleCIMessage );
  
  char handleDIAGMessage[] = "PsociGAsubspace: Invalid g_diag mat";
  local_g_diag->checkHandle( handleDIAGMessage );
  
  local_maxroots = local_vector_set->CurrentNumSoughtRoots(); // Actual number of root desired
  local_maxsubspace = local_vector_set->CurrentMaxRoots();    // should be SoughtRoots + maxubspace
  local_nsub = local_vector_set->CurrentExistingRoots();      // This is constant since we must preallocate 
  local_nbasis = local_vector_set->CurrentNumBasis();         // Should be maxsef

  // Compare to g_Hv assume g_diag is okay
  
  int type;
  int ndims;
  int dims[2];
  local_g_Hv->inquire( &type, &ndims, dims);
  if ( dims[0] != local_nbasis ) {
    cerr << "Nonconforming matrices in PsociGAsubspace: Hv" << dims[0] << "versus " << local_nbasis << endl;
    GA::error("Nonconforming matrices in PsociGAsubspace: Hv", dims[0]-local_nbasis ); 
  }
  
  local_g_diag->inquire( &type, &ndims, dims);
  if ( dims[0] != local_nbasis ) {
    cerr << "Nonconforming matrices in PsociGAsubspace:diag " << dims[0] << "versus " << local_nbasis << endl;
    GA::error("Nonconforming matrices in PsociGAsubspace:diag ", dims[0]-local_nbasis ); 
  }
  g_rnorms->inquire( &type, &ndims, dims);
  if ( dims[0] != local_nbasis  ) {
    cerr << "Nonconforming matrices in PsociGAsubspace:diag " << dims[0] << " " << dims[1] << "versus " << local_nbasis <<" local_maxroots" << endl;
    GA::error("Nonconforming matrices in PsociGAsubspace:diag ", dims[0]-local_nbasis ); 
  }
  
// Allocate arrays for the subspaces

  const int DATA_TYPE = C_DBL;
  const int MATRIX_NDIM = 2;
  const int VECTOR_NDIM = 1;
  dims[0] = local_maxsubspace;
  dims[1] = local_maxsubspace;
  
  if ( GA::nodeid() == 0 ) cout << "Create a PROJECTED GA of dims " << dims[0] << " by " << dims[1] << endl;
  
  const int  DISTRIBUTE_EVENLY = -1;
  int chunk[MATRIX_NDIM ];
  //chunk[0] = DISTRIBUTE_EVENLY; //TODO optimize these
  //chunk[1] = DISTRIBUTE_EVENLY;
  
  chunk[0] = dims[0]; //TODO optimize these
  chunk[1] = DISTRIBUTE_EVENLY; // Important: see diagSubsdpaceGOP

  g_s = GA::createGA( DATA_TYPE, MATRIX_NDIM, dims, (char *)"g_s overlap", chunk);
  g_vHv = GA::createGA( DATA_TYPE, MATRIX_NDIM, dims, (char *)"g_vHv ", chunk);

  g_evecs = GA::createGA( DATA_TYPE, MATRIX_NDIM, dims, (char *)"g_evecs subspace", chunk);


  chunk[0] =DISTRIBUTE_EVENLY;

  g_evals = GA::createGA( DATA_TYPE, VECTOR_NDIM, dims, (char *)"subspace evals", chunk);

  g_vHv->distribution( GA::nodeid(), local_vHvlo, local_vHvhi ); // Used by GOP methods
  
  if ( g_vHv->compareDistr( g_evecs ) ) {
    GA::error(" g_vHv not distributed as g_evecs ",-1);
  }
  
  g_s->zero();
  g_vHv->zero();
  g_evecs->zero();
  g_evals->zero();

#ifdef OLDGARNORM
  // One doesn't want to use these unless absolutely necessary
  dims[0]=local_nbasis;
  dims[1]=1;
  g_temp = GA::createGA( DATA_TYPE, MATRIX_NDIM, dims, (char *)"subspace temp ", NULL);
  g_temp2 = GA::createGA( DATA_TYPE, MATRIX_NDIM, dims, (char *)"subspace temp2 ", NULL);
  g_temp->zero();
  g_temp2->zero();
#endif

// Save time in incrementalSubspace and computernorms

  local_g_v->distribution( GA::nodeid(), vector_lo, vector_hi );
  local_g_Hv->distribution( GA::nodeid(), hvector_lo, hvector_hi );

  if ( vector_lo[0] != hvector_lo[0] || vector_lo[1] != hvector_lo[1] ||
	vector_hi[0] != hvector_hi[0] || vector_hi[1] != hvector_hi[1] ) {
        local_g_v->printDistribution();
        local_g_Hv->printDistribution();
        GA::error(" vector and Hvector not distributed the same ",-1);
   }

}

// In case the app forgets.......
/*
PsociGAsubspace::~PsociGAsubspace()
{
  destroySubspace(); // If previously destroyed then no harm is caused.
}
*/


GA::GlobalArray * PsociGAsubspace::fetchEvecsHandle()
{
   return ( g_evecs );
}

int PsociGAsubspace::numRoots() // subspace roots not necessarily SoughtRoots
{
  return( local_vector_set->CurrentExistingRoots() ); // Subspace vectors are always Existing_roots
} 

// Clearify what this is for ?
int PsociGAsubspace::numVectorBasis()
{
  return( local_nbasis ); // maxsef
}

//Returns the number of Rows since that will not be sparse. Columns are potentially shortened
int PsociGAsubspace::numCIBasis()
{
  return( local_nbasis );
}

// Destroy the current arrays created for the current subspace
void PsociGAsubspace::destroySubspace()
{
  if (GA::nodeid() == 0 ) cout << "Destroying PROJECTED GA subspace " << endl;
  g_s->destroy();
  g_vHv->destroy();
  g_evecs->destroy();
  g_evals->destroy();
#ifdef OLDGARNORM
  g_temp->destroy();
  g_temp2->destroy(); 
#endif
}

//Compute v(*)v and v(*)Hv - Hv is precalculated for all (existing set) roots
//Actual H computation not included yet.

/** The intent here is to preallocate all the subspace GA memory for the duration of the calculation.
    This is simple to do. the full space rank is maxsef and never changes. The subspace rank is EstimatedRoots
    but never gets any bigger than SoughtRoots+maxsubspace. So simply allocated the maximum possible size
    of the root index
*/

// GA-centric approaches at this stage are probably effective.
// TODO do not update entyire subspace everytime it is inefficient
// Collective
void PsociGAsubspace::generateSubspace()
{
  local_nsub = local_vector_set->CurrentExistingRoots(); //sought plus subspace
  
  /* Must use the PATCHED version of dgemm via matmul - yuk
     We could use dgemm BUT we would be doing lots of potentially extra work creating/destroying GA arrays
     We get these weird C-Fortran order issues creeping in.
  */

  int ailo,aihi,ajlo,ajhi;
  int bilo,bihi,bjlo,bjhi;
  int cilo,cihi,cjlo,cjhi;
  
  int n = local_nsub;
  int m = local_nbasis;
  
  // Strangeness in the settingof these values is related to weirdness in what GA does underneath
  // Check the hpctool google group list. Ha!

  ailo=0;
  aihi=n-1;
  ajlo=0;
  ajhi=m-1;
  
  bilo=0;
  bihi=m-1;
  bjlo=0;
  bjhi=n-1;
  
  cilo=0;
  cihi=n-1;
  cjlo=0;
  cjhi=n-1;
  
  double dOne=1.0, dZero=0.0;
  
  g_s->matmulPatch('N', 'T', &dOne, &dZero, local_g_v, ailo,aihi,ajlo,ajhi,
		   local_g_v, bilo,bihi,bjlo,bjhi, cilo,cihi,cjlo,cjhi );
  
  g_vHv->matmulPatch('N', 'T', &dOne, &dZero, local_g_v, ailo,aihi,ajlo,ajhi,
		     local_g_Hv, bilo,bihi,bjlo,bjhi, cilo,cihi,cjlo,cjhi );
}


// GA-centric approaches at this stage are probably effective.
/* In this case we simply want to do a partial update to the vHv/vSv matrices.  Here we
   specify in the argument list the NUMBER of new vectors in V and HV. These have been 
   orthonormalized and H*v already. The number of "existing_cols has bveen updated to reflect 
   the final set. With the value of NUM, we go back num entries and vHv/vVs those entries.

*/
//Collective


void PsociGAsubspace::generateIncrementalSubspace(int num)
{
  double dZero = 0.0;
  local_nsub = local_vector_set->CurrentExistingRoots(); //sought plus subspace
  int ld;
//  GA::GlobalArray * local_vec = local_g_v;

  /* SINCE the "new" vectors are already in the local_vector_set space, the value for
   existing_cols already include them. BUT, the incremental vHv space hasn;t yet been 
   updated. So, beginning at existing_cols, decrement "num" times and start from there.
   
   Not we currently only update the COLUMN and not the ROW
   I am not sure if we should explicitely symetrize the ROW
  */
  
  GA::sync(); // Important ! - the solution process can result in cores going rogue.

  int start = local_nsub - num;
  int end = local_nsub;
  
  // both are maxsef X roots. distributon of maxsef must be the same

/*
  int vilo[2], vihi[2];
  local_vec->distribution( GA::nodeid(), vilo, vihi );
  int vlow = vilo[0];
  int vhi  = vihi[0];
*/
  int vlow = vector_lo[0]; // pregenerated
  int vhi = vector_hi[0];

  int width = vhi - vlow + 1;
  
  int lo[2],hi[2];
  
  vector<double> veci( width, dZero ); 
  vector<double> hveci( width, dZero );
  vector<double> vecj( width, dZero ); // treat as scratch for vecj and Hvecj
  //vector<double> hvecj( width, dZero );

  vector<double> subv( local_nsub, dZero );
  vector<double> hsubv( local_nsub, dZero );
  
  int numElems = local_nsub;
  char op[] ="+";

//TODO check is a potential problem can occur if lo = -1

  for(int i=start; i< end; ++i ) {
    lo[0] = vlow; 
    hi[0] = vhi;
    lo[1] = i;
    hi[1] = i;
    ld = 1;
    
    local_g_v->get( lo, hi, &veci[0], &ld );
    local_g_Hv->get( lo, hi, &hveci[0], &ld );
    
    for ( int j=0; j< local_nsub; ++j ) 
      {
	lo[1] = j;
	hi[1] = j;
	local_g_v->get( lo, hi, &vecj[0], &ld );
        subv[j] = ddotPsoci( width, veci, ld, vecj, ld );

	local_g_Hv->get( lo, hi, &vecj[0], &ld );
        hsubv[j] = ddotPsoci( width, veci, ld, vecj, ld );
      }

      GA::gop( (double *) &subv[0], numElems, op );
      GA::gop( (double *) &hsubv[0], numElems, op );

/* TODO here we ONLY compute the unique new update. THis may be okay since the 
   diagonalization only grabs the triangle. But this could be misleading
*/
    
    lo[0] = 0;
    hi[0] = local_nsub - 1; 
    lo[1] = i;
    hi[1] = i;
    //ld = 1;

// Need to do this data-parallel
    if ( GA::nodeid() == 0 ) {
      g_s->put(lo, hi, &subv[0], &ld );
      g_vHv->put(lo, hi, &hsubv[0], &ld );
    }  
  }
  GA::sync();
}
/* GA doesn't provide a "patch" diag. As a result we cannot use the full
   maxsub * mabsub  g_vHv because that matrix is undiagonalizable. So we need to do the
   below patch copies
*/
//Collective
void PsociGAsubspace::diagSubspaceGA()
{
  GA::sync();
  int g_rank = GA::nodeid();
  
  local_nsub = local_vector_set->CurrentExistingRoots();
  
  int dims[ 2 ];
  dims[0] = local_nsub;
  dims[1] = local_nsub;
  
  GA::GlobalArray * g_patch_vHv = GA::createGA( C_DBL, 2, dims, (char *)"g_vHv patch subarray ", NULL);
  GA::GlobalArray * g_patch_v = GA::createGA( C_DBL, 2, dims, (char *)"g_v subarray ", NULL);
  GA::GlobalArray * g_patch_s = GA::createGA( C_DBL, 2, dims, (char *)"g_v subarray ", NULL);

  int lo[2], hi[2]; //vHv full size
  lo[0] = 0;
  hi[0] = local_nsub - 1;
  lo[1] = 0;
  hi[1] = local_nsub - 1;

  g_patch_vHv->copyPatch('N', g_vHv, lo, hi, lo, hi );
  g_patch_s->copyPatch('N', g_s, lo, hi, lo, hi );
  
  vector<double> evals( local_nsub, 0.0);
  g_evals->zero();
  
  g_patch_vHv->diagSeq( g_patch_s, g_patch_v, &evals[0] );
  g_evecs->copyPatch('N', g_patch_v, lo, hi, lo, hi );
  
#if 0
    vector<double>::iterator it;
    if ( GA::nodeid() == 0 ) { 
      for (it = evals.begin(); it != evals.end(); ++it) {
	cout << "DIAG: (local) ALL New Eigenvalues (sought+subspace) are " << (*it) << endl;
      }
    }
    cout << endl;
#endif   

    //push evals into a GA for subsequent rnorm computation
    if ( g_rank == 0 ) {
      int vlo = 0;
      int vhi = local_nsub - 1;
      int ldv = local_nsub;
      g_evals->put( &vlo, &vhi, &evals[0], &ldv );
    }  

    brdcstEvals( evals );

    g_patch_vHv->destroy();
    g_patch_v->destroy();
    g_patch_s->destroy();
}

//Collective but only root 0 does the print

void printCurrentEvals()
{
  if ( GA::nodeid() != 0 ) return;
/*
  int type;
  int dims[ 2 ];
  int idim;
  g_evecs->inquire( &type, &idim , dims );
  
  if ( num < 1 || num > dim[0] || num > dim[1] ) {
    GA::error("printCurrentEvals: erroneous value for num ",num);
  } 
  
  int lo[2],hi[2];
  int ld=1;
  
  lo[1] = 0;
  hi[1] = local_nsub - 1; // These are fixed and within class scope
  
  vector<double> evec( local_nsub, 0.0 );
  
  for(int i=0; i< num; ++i ) {
    lo[0] = iroot;
    hi[0] = iroot;
    lo[1] = 0;
    hi[1] = local_nsub - 1;
    ld = local_nsub;
    g_evecs->get( lo, hi, &evecs[0], &ld );
  }
*/
}

  // Collective - must run generate subspace and diag first 
  // On return send back all rnorms for the current set
  // Currently not really setup for parallel yet

  //GA++ doesn't seem to really suppoort my needs for this very well.

int PsociGAsubspace::computeRnormsGA( vector<double> & evals, vector<double> & rnorms )
{
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAsubspace: Invalid rnorm";
  g_rnorms->checkHandle( handleMessage );
#endif 
  
#ifndef OLDGASUBSPACE
  cerr << "computeRnormsGA can only be used by compiling the code with ";
  cerr << " -DOLDGASUBSPACE. This is not recommended though for perfomance reasons" << endl;
  GA::error("Must recompile the code to use computeRnormsGA",-1);
#endif
  g_temp->zero();
  g_temp2->zero();
  g_rnorms->zero();
  
  int numSoughtRoots = local_vector_set->CurrentNumSoughtRoots();
  
  int ld;
  double dOne=1.0;
  double dZero=0.0;
  
  /* from matmul.c
     m = *aihi - *ailo +1;
     n = *bjhi - *bjlo +1;
     k = *ajhi - *ajlo +1;
  */
  
  int lo[2],hi[2];
  int tlo[2],thi[2];
  
 // vector<double> evals( local_vector_set->CurrentExistingRoots(), dZero );
 // cout << "rnormsGA " << local_vector_set->CurrentExistingRoots() << " " << evals.size() << endl;
 // brdcstEvals( evals ); //Local to all

  for(int iroot=0; iroot< numSoughtRoots; ++iroot )
    {
      //get evecs for this root - replicated to all
      
      vector<double> evecs( local_nsub, dZero );
      lo[0] = iroot; //GADIAG forces this flipped order
      hi[0] = iroot; 
      lo[1] = 0;
      hi[1] = local_nsub - 1;
      ld = local_nsub; // applies to the local buffer
      if ( GA::nodeid() == 0 ) g_evecs->get( lo, hi, &evecs[0], &ld );
      GA::brdcst( &evecs[0], local_nsub*sizeof(double), 0) ;
      
      for(int isub=0; isub < local_nsub; ++isub )
	{
	  //g_temp->zero();
	  
	  lo[0] = 0;
	  hi[0] = local_nbasis - 1;
	  lo[1] = isub;
	  hi[1] = isub;
	  
	  tlo[0] = 0;
	  thi[0] = local_nbasis - 1;
	  tlo[1] = 0;
	  thi[1] = 0;
	  
       
	  g_temp->copyPatch('N',local_g_Hv, lo,hi,tlo,thi ); 
	  g_temp->scale( &evecs[ isub ] ); 
	  
	  g_temp2->addPatch(&dOne, g_temp2, tlo, thi, &dOne, g_temp, tlo, thi, tlo, thi );
	  
	  double scale = (-1.0)*evals[iroot]*evecs[isub];

	  g_temp->copyPatch( 'N', local_g_v, lo,hi,tlo,thi );
	  g_temp->scale( &scale );
	  
	  g_temp2->addPatch(&dOne, g_temp2, tlo, thi, &dOne, g_temp, tlo, thi, tlo, thi );
	}
      
       lo[0] = 0;
       hi[0] = local_nbasis - 1;
       lo[1] = iroot;
       hi[1] = iroot;
      
      g_rnorms->copyPatch( 'N', g_temp2, tlo, thi, lo, hi ); 
    }

  int rnum = fetchEigenvectorRnorms( rnorms );
  return( rnorms.size() );
}

/* A new method to take the current evecs and expand them into the full space.
   Reuse g_rnorms for this. The intent is to pass these to PsociGAbasis and replace
   the current set of existing roots with this new set of sought-size roots.
   This is a way to compress the subspace. NOTE: new H*v will be required so
   reset the appropriate parameters.
*/
// On output g_norms is overwritten

int PsociGAsubspace::expandEvecsGA()
{
  GA::error(" PsociGAsubspace::expandEvecsGA has been disabled use PsociGABasis methods ",-1);
/*
  char * handleMessage = "PsociGAsubspace: invalid scratch space";
  g_rnorms->checkHandle( handleMessage );

  int m = local_nbasis;
  int n = local_nsub;
  int k = local_nsub;
  double dOne=1.0;
  double dZero=0.0;
  double dNegate = -1.0;

  //cout << " M N K are " << m << " " << n << " " << k << endl;

  int ailo,aihi,ajlo,ajhi;
  int bilo,bihi,bjlo,bjhi;
  int cilo,cihi,cjlo,cjhi;
  
  ailo=0;
  aihi=n-1;
  ajlo=0;
  ajhi=m-1;
  
  bilo=0;
  bihi=n-1;
  bjlo=0;
  bjhi=n-1;
  
  cilo=0;
  cihi=n-1;
  cjlo=0;
  cjhi=m-1;

  // This will always grab all terms I think.....it is not a patch aproach but it shoud still work
  // g_rnorms->dgemm('N', 'N', m ,n, k, dOne, local_g_Hv, g_evecs, dZero );
  
  //g_rnorms->zero(); // May not be needed

  g_rnorms->matmulPatch('N', 'N', &dOne, &dZero, g_evecs, bilo,bihi,bjlo,bjhi,
			local_g_v, ailo,aihi,ajlo,ajhi, cilo,cihi,cjlo,cjhi ); 
  
*/
  return( local_nsub );
}

int PsociGAsubspace::expandEvecs( int num )
{
  GA::error(" PsociGAsubspace::expandEvecs has been disabled use PsociGABasis methods ",-1);

/*
  char * handleMessage = "PsociGAsubspace: invalid scratch space";
  g_rnorms->checkHandle( handleMessage );
  `
  if ( num > local_nsub ) {
     GA::error(" Nonesense value for num at expand ",num );
  }

  int m = local_nbasis;
  int n = num;
  int k = num;

  double dOne=1.0;
  double dZero=0.0;
  double dNegate = -1.0;

  int hvilo[2], hvihi[2];
  local_g_Hv->distribution( GA::nodeid(), hvilo, hvihi );
  int hvlow = hvilo[0];
  int hvhi  = hvihi[0];

  int ailo,aihi,ajlo,ajhi;
  int bilo,bihi,bjlo,bjhi;
  int cilo,cihi,cjlo,cjhi;
  
  ailo=0;
  aihi=n-1;
  ajlo=0;
  ajhi=m-1;
  
  bilo=0;
  bihi=n-1;
  bjlo=0;
  bjhi=n-1;
  
  cilo=0;
  cihi=n-1;
  cjlo=0;
  cjhi=m-1;

  cout << "M and N are " << n << " " << m << endl;
  g_rnorms->matmulPatch('N', 'N', &dOne, &dZero, g_evecs, bilo,bihi,bjlo,bjhi,
			local_g_v, ailo,aihi,ajlo,ajhi, cilo,cihi,cjlo,cjhi ); 
*/
  return( local_nsub );
}

//Collective call but only root=0 does the get
int  PsociGAsubspace::brdcstEvals( vector<double> & evals )
{
  local_nsub = local_vector_set->CurrentExistingRoots();
  int vlo = 0;
  int vhi = local_nsub - 1;
  int ldv = local_nsub;
  
#ifdef DETAILEDCHECK
  if ( evals.size() < local_nsub ) {
    cerr << "Bad evals array size " << evals.size() << endl;
    GA::error("Bad evals array size ", evals.size() );
  }
#endif

  GA::sync(); // Important
  if ( GA::nodeid() == 0 ) g_evals->get( &vlo, &vhi, &evals[0], &ldv );

  GA::brdcst( &evals[0], local_nsub*sizeof(double), 0) ;

  //cout << "extra sync after brdcst " << endl;
  //GA::sync();
  return( local_nsub );
}


// The sum of SoughtRoots + maxsubspace
int PsociGAsubspace::fetchMaxSubspace()
{
  return( local_maxsubspace );
}

int PsociGAsubspace::fetchExistingSubspace()
{
  return( local_vector_set->CurrentExistingRoots() );
}

/** Compute r^t*r for all current existing vectors
    Let the calling app decide which ones are important
*/
//COLLECTIVE
int PsociGAsubspace::fetchEigenvectorRnorms( vector<double> & rnorms )
{
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAsubspace:getrnorms";
  g_rnorms->checkHandle( handleMessage );
#endif
  
  int existing_cols = local_vector_set->CurrentExistingRoots();
  int local_soughtroots = local_maxroots;
  int global_rows = local_vector_set->CurrentNumBasis();

#ifdef DETAILEDCHECK
    if ( rnorms.size() > 0 ) {
    cout << "getAllresidual norms: Warning: expected an empty norms array. Emptying and proceeding " << rnorms.size() << endl;
    }
  if ( rnorms.size() !=  local_soughtroots ) {
    GA::error(" Bad rnorms size " , rnorms.size() );
  }
#endif

  //rnorms.resize( local_soughtroots, 0.0 );
  
  // Will place temporary values in the local array double V[roots] then global sum the whole thing
  // Or we could stick results into a GA vector of length nroots ?
  // Only process existing_rows. If/when new vectors are added to the groups, then update existing to reflect this.
  
  // Note distribution provides dims in C-indexing style (0,n-1)
  
  int lo[2], hi[2];

  int rootlo = 0;
  int roothi = local_soughtroots-1;
  
  int basislo = 0;
  int basishi = global_rows-1;
  int normstatus = 0;
  
  double precision = local_vector_set->fetchMinPrecision();
  double value;
  
  lo[0] = basislo;
  hi[0] = basishi;
  
  for (int i=rootlo; i <= roothi; ++i )
    {
      lo[1]=i; // Process residual at a time
      hi[1]=i;
      value = sqrt(g_rnorms->ddotPatch('N', lo, hi, g_rnorms, 'N', lo, hi ));
      if ( value < precision ) { // low Frobenius norm okay for these - in fact trending->zero is a good sign
        normstatus = 1;
      }
      rnorms[i] = value;
    } 
  return( normstatus );
}

void PsociGAsubspace::printVectorContents( vector<double> & rnorms )
{
  if ( GA::nodeid() != 0 ) return;
  cout << "Printing out vector (rnorm) contents " << endl;
  vector<double>::iterator it;
  for(it = rnorms.begin(); it != rnorms.end(); ++it ) {
    cout << (*it) << " " ;
  }
  cout << endl;
}
  // Collective - must run generate subspace and diag first 
  // On return send back all rnorms for the current set
  // Currently not really setup for parallel yet
  // This version uses local BLAS for efficiency

/* The performance and memory footprint of the GA method
   were far to much for the subspace operations. As such, 
   it was neccesssary to go back to the original (cidbg4) model
   which was a data-parallel approach using local BLAS
*/

int PsociGAsubspace::computeRnorms( vector<double> & evals, vector<double> & rnorms )
{
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAsubspace: Invalid rnorm";
  g_rnorms->checkHandle( handleMessage );
#endif
  
  //g_rnorms->zero(); // May not be needed
  
  int numSoughtRoots = local_vector_set->CurrentNumSoughtRoots();
  
  int ld;
  double dZero=0.0;
  double dNegate = -1.0;

  /* from matmul.c
     m = *aihi - *ailo +1;
     n = *bjhi - *bjlo +1;
     k = *ajhi - *ajlo +1;
  */

  int lo[2],hi[2];
  int width;
  double value;
  
  //vector<double> evals( local_nsub, dZero );
  //vector<double> evals( local_vector_set->CurrentExistingRoots(), dZero );

//cout << "rnorms regulat " << local_vector_set->CurrentExistingRoots() << " " << evals.size() << endl;
 // brdcstEvals( evals ); //Local to all
  
  // both are maxsef X roots. distributon of maxsef must be the same
/*
  int vilo[2], vihi[2];
  local_vector_set->fetchHandle()->distribution( GA::nodeid(), vilo, vihi );
  int vlow = vilo[0];
  int vhi  = vihi[0];
*/
  int vlow = vector_lo[0];
  int vhi = vector_hi[0];
  
  int hvlow = vector_lo[0];
  int hvhi  = vector_hi[0];

  lo[0] = vlow;
  hi[0] = vhi;
  width = hvhi - hvlow +1;
  for(int iroot=0; iroot< numSoughtRoots; ++iroot )
    { 
      lo[1] = iroot;
      hi[1] = iroot;
      
      //get evecs for this root - replicated to all
      vector<double> evecs( local_nsub, dZero );
      lo[0] = iroot; //GADIAG forces this flipped order
      hi[0] = iroot;
      lo[1] = 0;
      hi[1] = local_nsub - 1;
      ld = local_nsub;
      if ( GA::nodeid() == 0 ) g_evecs->get( lo, hi, &evecs[0], &ld );
      GA::brdcst( &evecs[0], local_nsub*sizeof(double), 0) ;
      
      vector<double> tvec( width, dZero );
      vector<double> hvec( width, dZero);
      
      lo[0] = hvlow;
      hi[0] = hvhi;
      double evalue = evals[iroot];
      
      for(int isub=0; isub < local_nsub; ++isub )
        {
	  lo[1] = isub;
	  hi[1] = isub;
	  ld = 1;
          double evector = evecs[ isub ];
	  
          local_g_Hv->get( lo, hi, &tvec[0], &ld );  

          daxpyPsoci( width,  evector, tvec, ld, hvec, ld);

          local_vector_set->fetchHandle()->get( lo, hi, &tvec[0], &ld );
          value = dNegate *  evector * evalue; 

          daxpyPsoci( width, value, tvec, ld, hvec, ld);
        } 
      lo[1] = iroot;
      hi[1] = iroot;
      g_rnorms->put( lo, hi, &hvec[0], &ld );
    }
  int rnum = fetchEigenvectorRnorms( rnorms );
  return( rnorms.size() );
}

/* Return the ||rnorm|| for the root; iroot and a return value of -1 it not converged, else >=0
 */

// iroot in comes as Fortran order (1,n)
// Compute a single rnorm based on the imput hint

int PsociGAsubspace::computeRnormSingle( int iroot, vector<double> & evals, double & rnormvec )
{
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAsubspace: Invalid rnorm";
  g_rnorms->checkHandle( handleMessage );
#endif
  
  int numSoughtRoots = local_vector_set->CurrentNumSoughtRoots();
  
//#ifdef DETAILEDCHECK
  if ( iroot > numSoughtRoots ) {
     cout << "computeRnormSingle: iroot > numSoughtRoots "<<iroot<<endl;
     GA::error("computeRnormSingle: iroot > numSoughtRoots ", iroot );
  }
//#endif
  
  int ld;
  const double dZero=0.0;
  const double dNegate = -1.0;
  
  rnormvec = dZero;

  int liroot = iroot - 1;
  int lo[2],hi[2];
  int width;
  double value;
  
  int normstatus = 0;

/* This entry was changed because when using the "single" version we desire
   that at the end of a iteration cycle to fully calculate all the final rnorms.
   However, the value of local_nsub got incremented in diagSubspace which would occur
   AFTER this step. Thus change to the more correct "CurrentExistingRoots()" call
*/
 // Deleted old code that incorrectly openened up new evals space
 // vector<double> evals( local_vector_set->CurrentExistingRoots(), dZero );
 // cout << "rnorm Single " << local_vector_set->CurrentExistingRoots() << " " << evals.size() << endl;
 // brdcstEvals( evals ); //Local to all - hides the evals in driver
/*
  int vilo[2], vihi[2];
  local_vector_set->fetchHandle()->distribution( GA::nodeid(), vilo, vihi );
  int vlow = vilo[0];
  int vhi  = vihi[0];
*/

  int vlow = vector_lo[0];
  int vhi = vector_hi[0];

  lo[0] = vlow;
  hi[0] = vhi;
  lo[1] = liroot;
  hi[1] = liroot;
  width = vhi - vlow + 1;
  
  vector<double> evecs( local_nsub, dZero );
  lo[0] = liroot; //GADIAG forces this flipped order
  hi[0] = liroot;
  lo[1] = 0;
  hi[1] = local_nsub - 1;
  ld = local_nsub;

// we could consider brdcasting evas as well
  if ( GA::nodeid() == 0 ) g_evecs->get( lo, hi, &evecs[0], &ld );
  GA::brdcst( &evecs[0], local_nsub*sizeof(double), 0) ;
  
  vector<double> tvec( width, dZero );
  vector<double> hvec( width, dZero );
  double evalue = evals[liroot];

  lo[0] = vlow;
  hi[0] = vhi;
  ld = 1;
//  if ( vlow != < 0 ) {  This is unlikely given it's for a fullspace coordinate
  for(int isub=0; isub < local_nsub; ++isub )
    {
      lo[1] = isub;
      hi[1] = isub;
      double evector = evecs[ isub ];
      
      local_g_Hv->get( lo, hi, &tvec[0], &ld );  
      
      daxpyPsoci( width, evector, tvec, ld, hvec, ld);

      local_vector_set->fetchHandle()->get( lo, hi, &tvec[0], &ld );
      value = dNegate *  evector * evalue; 
      
      daxpyPsoci( width, value, tvec, ld, hvec, ld);
    }
  
    rnormvec = ddotPsoci( width, hvec, ld, hvec, ld );
//    } // actually had local data

//assemble distributed rnorm

  char op[] = "+";
  int numElems  = 1;
  GA::gop( (double *) &rnormvec, numElems, op );

  rnormvec = sqrt( rnormvec );

//  double precision = local_vector_set->fetchMinPrecision();
//  if ( rnormvec < precision ) { // low Frobenius norm okay for these - in fact trending->zero is a good sign

  lo[0] = vlow;
  hi[0] = vhi;
  lo[1] = liroot;
  hi[1] = liroot;
  ld = 1;
//  if ( vlow != < 0 ) {
  g_rnorms->put( lo, hi, &hvec[0], &ld ); //put it in more often than not this will be correct 
//  }
  return( normstatus );
}

/* GA doesn't provide a "patch" diag. As a result we cannot use the full
   axsub * mabsub  g_vHv because that matrix is undiagonalizable. So we need to do the
   below patch copies
*/

/* The GA diag seems to not work too well. Since it does al of its work on a single node, 
   It seems better to do it myself, not but the price of al the etxra A calls and use a LAPACK
   routine for a local diagonalization. 
   
   This will soon be a bottle neck but for now this should be okay
*/
//Collective
void PsociGAsubspace::diagSubspace()
{
  int g_rank = GA::nodeid();
  

  local_nsub = local_vector_set->CurrentExistingRoots();

  int lo[2], hi[2]; //vHv full size
  lo[0] = 0;
  hi[0] = local_nsub - 1;
  lo[1] = 0;
  hi[1] = local_nsub - 1;
  
/* This is a simple procedure. The reduced matrices are already constructed.
     Simply read into local memory on one core and call SSYGV
*/

  int itype = 1; // A*x = B*x*lambda
  char jobz = 'V';
  char uplo = 'L'; //doesn;t really matter since vHv is currently symmetric and dense
  int lwork = 50 * local_nsub; //TODO Need to optimize this
  int info;

/* Rework this. For big problems the diag time is dominated by GA. So instead of
   node 0 simply fetching every entry, do a fetch data-parallel, then use optimized
   gops to assemble the data
*/

  vector<double > local_s(local_nsub *  local_nsub, 0.0);
  vector<double > local_vHv(local_nsub * local_nsub, 0.0);
  vector<double > work(lwork, 0.0); // could be optimized abit. 
  vector<double> evals(local_nsub, 0.0);
  
  int ld;
  double tempTime;
  
  if ( GA::nodeid() == 0 )
    {
      ld = local_nsub;
      tempTime = psociTime();
      g_s->get( lo, hi, &local_s[0], &ld );
      
      ld = local_nsub;
      g_vHv->get( lo, hi, &local_vHv[0], &ld );
      
      // SSYGV is the LAPACK equivelent to RSG
      // NOTE we must worry about packing order and compress into triangular shape
      // ssygv | dsygv (iopt, a, lda, b, ldb, w, z, ldz, n, aux, naux);
      // cout << "DIAG GET TIME is " << psociTime() - tempTime << endl;
      
      dsygvPsoci( itype, jobz, uplo, local_nsub,  local_vHv, local_nsub, local_s, local_nsub, evals, work, lwork, info );
      
#if 0
      for(int i=0; i< local_nsub; ++i ) {
	cout << "EVALS " << evals[i] << endl;
      }
      vector<double>::iterator it;
      if ( GA::nodeid() == 0 ) { 
	for (it = evals.begin(); it != evals.end(); ++it) {
	  cout << "DIAG: (local) ALL New Eigenvalues (sought+subspace) are " << (*it) << endl;
	}
      }
      cout << endl;
#endif
      
      //push evals into a GA for subsequent rnorm computation
      int vlo = 0;
      int vhi = local_nsub - 1;
      int ldv = local_nsub;
      g_evals->put( &vlo, &vhi, &evals[0], &ldv );
      
      //Need to push local eigenvectors into GA
      lo[0] = 0;
      hi[0] = local_nsub - 1;
      lo[1] = 0;
      hi[1] = local_nsub - 1;
      ld = local_nsub;
      g_evecs->put( lo, hi, &local_vHv[0], &ld );      
    } // GA::nodeid
  // brdcstEvals( evals );
}

//NEW VERSION Collective
// Substitutes GOPS for lots of remote GA gets

void PsociGAsubspace::diagSubspaceGOP()
{
  int g_rank = GA::nodeid();  
  char op[] = "+";
  
  local_nsub = local_vector_set->CurrentExistingRoots();
  
  int lo[2], hi[2]; //vHv full size
  /* old way
    lo[0] = 0;
    hi[0] = local_nsub - 1;
    lo[1] = 0;
    hi[1] = local_nsub - 1;
  */
  
  /* This is a simple procedure. The reduced matrices are already constructed.
     Simply read into local memory on one core and call DSYGV
  */
  
  int itype = 1; // A*x = B*x*lambda
  char jobz = 'V';
  char uplo = 'L'; //doesn;t really matter since vHv is currently symmetric and dense
  int lwork = 50 * local_nsub; //TODO Need to optimize this
  int info;
  double dZero = 0.0;
  
  vector<double > local_s(local_nsub *  local_nsub, dZero ) ;
  vector<double > local_vHv(local_nsub * local_nsub, dZero);
  vector<double > work(lwork, dZero); // could be optimized abit. 
  vector<double> evals(local_nsub, dZero) ;

  int index = PsociGAsubspace::localPatchFetch( lo, hi, local_nsub );


#if 0
  cout << local_nsub << " status " << index << " "<< lo[0]<<" "<<lo[1]<<" "<<hi[0]<<" "<<hi[1] << endl;
  g_s->printDistribution();
  g_vHv->printDistribution();
  g_evecs->printDistribution();
  g_evals->printDistribution();
#endif
  
  int ld = local_nsub;
  if ( index != -1 ) {
    g_s->get( lo, hi, &local_s[index], &ld );
    g_vHv->get( lo, hi, &local_vHv[index], &ld );
  }

// Leverage high speed collectives instead of GA gets/puts

  int numElems  = local_nsub * local_nsub;
  GA::gop( (double *) &local_s[0], numElems, op );
  GA::gop( (double *) &local_vHv[0], numElems, op );
  
  double tempTime = psociTime();
 
  dsygvPsoci( itype, jobz, uplo, local_nsub,  local_vHv, local_nsub, local_s, local_nsub, evals, work, lwork, info );
  
#if 0
  for(int i=0; i< local_nsub; ++i ) {
    cout << GA::nodeid() << " EVALS " << evals[i] << endl;
  }
  vector<double>::iterator it;
  if ( GA::nodeid() == 0 ) { 
    for (it = evals.begin(); it != evals.end(); ++it) {
      cout << "DIAG: (local) ALL New Eigenvalues (sought+subspace) are " << (*it) << endl;
    }
  }
  cout << endl;
#endif
 
/* This seems okay- I didn;t want to put a sync here but on occasion, the calling app 
   will pruintout the evals and sometimes not all nodes here  were finished.
   Since, in the end, the values are correct I leave this as it is.
*/

  if ( index != -1 ) {
    ld = local_nsub;
    g_evecs->put( lo, hi, &local_vHv[index], &ld );
    g_evals->put( &lo[1], &hi[1], &evals[index], &ld );
    }

  // brdcstEvals( evals ); 
}

int PsociGAsubspace::fetchExistingSubspaceSize()
{
  return( existing_subspace_size );
}

// A local version that ( for these cases ) works better than using MKL or ACML
inline void  PsociGAsubspace::daxpyLocal( int width, double value, double * tvec, double * hvec)
{
  for(int i=0; i< width; ++i ) {
    hvec[i] += tvec[i] * value;
  }
}

/* Inside of diagSubspaceGOP we attemp to save effort by keeping all GA fetch's local
   However, the current subspace size <= maximum subspace size. SO we need to detemrine
   who and how much should actually attemp to fetch and process data 
   
   for now we assume rows and columns behave similarly
*/

// Have this method return the index for the first insertion point ( or -1 )
inline int PsociGAsubspace::localPatchFetch( int * lo, int * hi, int maxvalue )
{
  int shiftmax = maxvalue - 1;
  int start_index = -1;
  
  int vlo1 = local_vHvlo[0];
  int vlo2 = local_vHvlo[1]; // Check both since distribution can be varied by the calling app
  int vhi1 = local_vHvhi[0];
  int vhi2 = local_vHvhi[1];
  
  if ( (vlo1 >= 0 && vlo1  <=  shiftmax) && (vlo2 >= 0 && vlo2 <= shiftmax) ) {
    //if ( vhi1 > shiftmax ) vhi1 = shiftmax;
    //if ( vhi2 > shiftmax ) vhi2 = shiftmax;
    lo[0] = vlo1;
    lo[1] = vlo2;
    hi[0] = min( shiftmax, vhi1);
    hi[1] = min( shiftmax, vhi2);
    start_index = vlo1 * maxvalue + vlo2; // helps the calling program 
  } else {
    lo[0] = -1;
    lo[1] = -1;
  }
  return( start_index );
}


