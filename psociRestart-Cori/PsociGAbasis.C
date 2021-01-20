/** TODO
  The usage of ddot versus sqrt(ddot) is cumbersome. Need to clean this up
*/

/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A rudimentary class managing the manipulation and processing of
              vectors in an orthonormal CI basis (or otherwise unit metric) stored in Global Arrays
              data are presumed stored as (basis, roots ) ( basis == rows, roots == cols );

              These are generally Collective calls. 

              These methods process sets of vectors (and H*v products). 
 History:

              Sept 2011. Preliminary validation completed

**************************************************************************************/
/**
 *   @file PsociGAbasis.C
 *
 */

#include "PsociGAbasis.hpp"
#include "PsociBlasLapack.hpp"

using namespace std;

/** This constructor assumes the number of desired cols == number of current cols AND that the
    GA has (at least) some representation of the vectors. They need not be orthonormal at this time.
    We can fix them later
*/

//Collective
PsociGAbasis::PsociGAbasis( int numSoughtRoots, PsociGAhamiltonian * hamiltonian, GA::GlobalArray * g_array,  GA::GlobalArray * g_Harray )
{
  GA::SERVICES.sync();
  //GA::sync();
  local_hamiltonian = hamiltonian; //useful for procond and H*v steps
  
  local_ga = g_array;   // the set of vectors
  local_Hga = g_Harray; // the set of related matrix vector products: precision ?
  
  char handleMessage[] = "PsociGAbasis-2";
  local_ga->checkHandle( handleMessage );
  
  char handleHMessage[] = "PsociGAbasis-2:Hv";
  local_Hga->checkHandle( handleHMessage );
  
  int dims[MAX_GA_VECTOR_DIMS];
  int ndim = local_ga->ndim();
  if ( ndim > MAX_GA_VECTOR_DIMS ) {
    cerr << "ndim is " << ndim << " and cannot be greated than " << MAX_GA_VECTOR_DIMS << endl;
    GA::Terminate();
  }

  local_ga->distribution( GA::nodeid(), vector_lo, vector_hi );
  local_Hga->distribution( GA::nodeid(), hvector_lo, hvector_hi );
  int type;
  local_ga->inquire( &type, &vector_ndim, vector_dims );

// No need to check conformace subspace code does that
  
  global_rows = vector_dims[0];
  global_cols = vector_dims[1]; // numSoughtRoots + maxsubspace

  existing_cols = numSoughtRoots; // eventually includes numSoughtRoots+present subspace(variable)
  actual_roots_sought = numSoughtRoots; // This remains FIXED unless manually resized.
  converged_cols = 0; // Okay at constructor always start from scratch for now
  desired_cols = COLS_DESIRED_UNSPECIFIED; //Triggered only when appending
  
  present_matrixProducts = 0;
  required_matrixProducts = existing_cols;
  
  if ( existing_cols > global_cols ) {
    cerr << "Existing rows > Global cols: Existing: " << existing_cols << " Global: " << global_cols << endl;
    cerr << "Though this is a recoverable error, I am aborting so you can double check what you are doing " << endl;
    GA::Terminate();
  }

/* Create a fresh g_array to act as a preconditioning scratch area
*/
  
  int DATA_TYPE = C_DBL; // 
  int DATA_NDIM = 2;
  dims[0] = global_rows;
  dims[1] = 1;
  int chunk[2];
 
  chunk[1] = 1;
  chunk[0]=-1;
  
  char local_title[] = "scratch vector array ";
  g_scratch = GA::createGA( DATA_TYPE, DATA_NDIM, dims, (char *) local_title, chunk);
  g_scratch->inquire( &type, &scratch_ndim, scratch_dims );

  l_scatter_products = false; 

}

void PsociGAbasis::setProducts( bool scatter_products )
{
   l_scatter_products = scatter_products;
}

bool PsociGAbasis::getProducts()
{
   return(l_scatter_products);
}

/*
PsociGAbasis::~PsociGAbasis()
{
  g_scratch->destroy(); 
}
*/


GA::GlobalArray * PsociGAbasis::fetchHandle(){
  return( local_ga );
}

GA::GlobalArray * PsociGAbasis::fetchHvHandle(){
  return( local_Hga );
}

// Use full if beginning subdspace work with preexisting H*v terms
void PsociGAbasis::setpresent_matrixProducts( int l_present_matrixProducts )
{
  present_matrixProducts = l_present_matrixProducts;
}

double PsociGAbasis::fetchMinPrecision()
{
  return( MIN_PRECISION );
}

//Local call
int PsociGAbasis::CurrentExistingRoots()
{
  return( existing_cols );
}

int PsociGAbasis::CurrentMaxRoots()
{
  return( global_cols );
}

int PsociGAbasis::CurrentNumSoughtRoots()
{
  return( actual_roots_sought );
}

int PsociGAbasis::CurrentNumBasis()
{
  return( global_rows );
}

int PsociGAbasis::CurrentReqHVProducts()
{
  return( required_matrixProducts );
}

int PsociGAbasis::CurrentExistingHVProducts()
{   
  return( present_matrixProducts);
} 


/* In this routine we adjust the size of the input GA to reflect the apps choice of
   only having a number of roots ==  desired_cols. When desired cols < local_cols
   we need to go ahead and simply create them. Note this is done (for now) by simply 
   adding random numbers to augmenting vectors. We will need to do a 
   modified gram-schmidt/normalization step afterwards.
   
   In any event we flush to ZERO all vector outside of desired_set iff < global_cols
*/

// Collective 
// timer wrapper 
// setCols is the TOTAL number of cols to exist in the GA
// We assume RESIZE is a user specified change and update existing_roots accordingly;

// This only changes the SUBSPACE size not the sought number of roots
int PsociGAbasis::resizeVectors( int setCols, pair<int,double> & time  )
{
  int rank = GA::nodeid();
  double timein = psociTime();
  int status = resizeVectors( setCols );
  time.first = rank;
  time.second = psociTime()-timein;
  return( status );
}

//Collective
/** The values for desired and existing rows have already been set.
    Either by the first constructor ( which sets existing == desired ) in which 
    case you cannot make any changes OR through the second constructor where
    sought_cols < actual_cols. Here one can add up to actual_col worth of extra
    vectors.
*/
/* On input specify the number of cols to use. On exit the GA is irreversibly changed
   To have fewer non-zero vectors - OR - new vectors are added/orthonormalized.
*/

// reset actual_roots_sought
// setCols is the TOTAL number of rows to exist in the GA
int PsociGAbasis::resizeVectors( int setCols )
{
  desired_cols = setCols;
  
#ifdef DETAILEDCHECK
//  if ( desired_cols > global_cols ) {
//  cerr << "Cannot add cols to this GA. It is already full: GA cols = " << global_cols;
//  cerr << " Requested resize is " << desired_cols << endl;
//  return( add_status );
//  }
//  if ( desired_cols == global_cols ) {
//  cout << "No cols changed for this GA. desired_cols = GA cols: GA cols  = " << global_cols;
//  cout << " Requested resize is: No action performed " << desired_cols << endl;
//  add_status = 0;
//  return( add_status );
//  }
#endif
  // Now we either shrink the GAs by zeroing out the highest numbered vectors OR augment
  // Cannot simply zero out higher contributions you lose too much in the optimization
  
  int basislo = 0;
  int basishi = global_rows-1;
  int lo[2], hi[2];
  //lo[1] = -1;
  //hi[1] = -1;
  lo[0] = basislo;
  hi[0] = basishi;
  
  if ( desired_cols < global_cols ) {
    
    // Two possibilities 1) Shrink GA starting at the back - resert actual number of roots . orthonormalization not performed;
    //if (GA::nodeid() == 0 ) cout << "COMPRESS VECTOR SPACE " << existing_cols <<" "<<desired_cols << " "  << global_cols << endl;
    
    if ( desired_cols < existing_cols ) {
      for( int i=(existing_cols -1); i > (desired_cols-1); --i )
	{
	  lo[1] = i; //Delete Ith col 
	  hi[1] = i;
	  //cout << "Eliminating(zeroing out) the col i " << (i+1) << endl; 
	  local_ga->zeroPatch( lo, hi );
	  --existing_cols;
          //should not change this actual_roots_sought = existing_cols;
	}
    }
    
    // 2) Expand existing set of GA vectors.
    
    //Now need to ADD vectors to the set. You should orthonormalize after this in the calling APP 
    //Fill vectors with random numbers for now - what else could one do ?
    /*
      For now use a simple RAND based approch. RAND is okay here since true
      randomness is not required..We only want some kind of initial guess
      Have every process fill its own local buf with random numbers then PUT back to the GA
    */
    
    if ( desired_cols >  existing_cols ) 
      {
	int me = GA::nodeid();
	int mylo[2], myhi[2];
/*
	local_ga->distribution( me, mylo, myhi );
*/
        mylo[0] = vector_lo[0];
        myhi[0] = vector_hi[0];
        mylo[1] = vector_lo[1];
        myhi[1] = vector_hi[1];
	
	int lowbas = mylo[0];
	int hibas = myhi[0];
	int lowrt = mylo[1];
	int hirt = myhi[1];
	
        
	int ld = 1;
	
	//Watch out for the damn index scheme switch
	for (int i = existing_cols; i < desired_cols; ++i) 
	  {
	    //skip the check on -1 for now
	    
	    if ( ((i-1) <= hirt) && ((i-1) >= lowrt) ) { //Then we own part of i
	      vector<double> buf;
	      fetchPseudoRandom( (hibas-lowbas+1), i*(GA::nodeid()+3), buf ); //Grab a set of random numbers
#if 0
	      vector<double>::iterator it;
	      for(it=buf.begin(); it!=buf.end(); ++it) {
		cout << "Random value is " <<  (*it) << endl;
	      }
#endif
	      if (buf.size() != (hibas-lowbas+1) ) {
		cerr << "Error: Unexpected size of Random number array: Aborting " << endl;
		return(-1);
	      }
	      
	      //Shove individual vector into the new vector slot.
	      mylo[1] = i;
	      myhi[1] = i;
	      local_ga->put(mylo,myhi,&buf[0],&ld);        
	    }
	  }
        existing_cols = desired_cols;
      }
    
    
    local_Hga->zero(); // Need to recalculate them all
#ifdef NEWORTHOALL
    int or_status = orthonormalizeAppendedVectors(0 , desired_cols);
#else
    int or_status = orthonormalizeAllVectors();
#endif
    
    // Now we need to compute new H*v for both conditions
    // Or since we orthogonaliuzed all vecs we recompute all H*v YES for now.
    // worry about redundant work later.
    
    present_matrixProducts = 0; //Recalculate all of them
    required_matrixProducts = existing_cols; // Restart subspace >= soughtNum
    
    if ( computeHvProducts( present_matrixProducts, required_matrixProducts ) < 0 ){
      cerr << "computeHvProducts failed " << endl;
      GA::error("computeHvProducts failed ",GA::nodeid() );
    } // Of course for now it will NEVER trigger a fail by construction 
    
    //present_matrixProducts = actual_roots_sought;
    present_matrixProducts = existing_cols;
    required_matrixProducts = 0;
    
    return( or_status );
  }
  cerr << "Should never fall out bottom " << endl;
  return( -1 );
}

// setRows is the TOTAL number of rows to exist in the GA
int PsociGAbasis::appendVectors( GA::GlobalArray * g_vectors, pair<int,double> & time  )
{
  int rank = GA::nodeid();
  double timein = psociTime();
  int status = appendVectors( g_vectors );
  time.first = rank;
  time.second = psociTime()-timein;
  return( status );
}

// A convenience method to handle appending the contents of the "internal" object g_scratch
int PsociGAbasis::appendScratchVector()
{
  int num = appendVectors( g_scratch );
  return( num );
}

/* A routine that takes a GA of size (1:nbasis, 1:nroots) and appends it to the current set of vectors g_v(nbasis,nroots)
   Currently no contraction is performed if no more space exists in g_v. Simply returns a failure for now.
*/

/* This method is used to compress the subspace back to size = num It does this by grabbing
   the current set of expanded evecs from PsociGAsubspace and simply selecting those as the new vectors.
   Of course, fresh H*v will be needed.  Vectors that are deemed CONVERGED should not be replaced but
   that wil be dealt with later. Then for each restart of the subspace we simply recompute the "converged"
   vectors which should converged in one step.
*/

/* Take NUM number of eigenvectors and expand into the full space then use those 
   as the complete set of subspace expansion vectors 
*/

//TODO some check that the number of requested Evectors expanded makes sense

int PsociGAbasis::expandAndReplaceVectors( GA::GlobalArray * g_evecs, int num )
{
#ifdef DETAILEDCHECK
  if( num > CurrentExistingRoots() ) {
    GA::error("expandAndReplaceVectors: num is wrong", num);
  }
#endif

  // The plan here is to expand these eigenvectors into the existing V and Hv space, 
  // local_ga = g_array;   // the set of vectors
  // local_Hga = g_Harray
  
  int local_nbasis = CurrentNumBasis();
  int local_nsub = CurrentExistingRoots();
  
  double dZero = 0.0;
  double value;
  
  /* With this new algorithm, local_nsub should never get more than 3-4*numSoughtRoots ~ 200
     So at most we might need to fetc 100 roots * 200 subspace ~ 2,000 * 8 = 16kB of data 
     This is per core and could be alot of GA latencies. So let's simply grab them ALL and brdcst
     which is probably optimized for all machines
  */
  
  int rld, rlo[2], rhi[2];
  vector<double> root( local_nsub * num );
  // in the format root[num][local_nsub];
  
  if ( GA::nodeid() == 0 ) {
    rlo[1]=0;
    rhi[1]=local_nsub-1;
    rlo[0]=0;
    rhi[0]=num-1;
    rld = local_nsub;
    g_evecs->get( rlo,rhi, &root[0], &rld );
  }
  GA::brdcst( &root[0], local_nsub*num*sizeof(double), 0);
  
  // DO THIS IN DATA-PARALLEL
  
  int vilo[2], vihi[2];
/*
  local_ga->distribution( GA::nodeid(), vilo, vihi );
*/
  int vlow = vector_lo[0];
  int vhi  = vector_hi[0];
  int vld;
  
  if ( vector_lo[0] < 0 || vector_lo[1] < 0 ) {
     cerr << "expand NO DATA ON ME " << GA::nodeid() << endl;
     GA::Terminate();
  }
  
  int width = vhi - vlow + 1;
  vector<double> vec( width, dZero ); //move these after initial validation
  vector<double> hvec( width, dZero );

  int index;

  vld = 1;
  vilo[0] = vlow;
  vihi[0] = vhi;
  for(int iroot=0; iroot<num; ++iroot )
    {
      for(int j=0; j< local_nsub; ++j ) 
	{
	  index = iroot * local_nsub + j; 
	  value = root[index]; //[root][local_nsub]
	  //vilo[0] = vlow;
	  //vihi[0] = vhi;
	  vilo[1] = j; //reuse
	  vihi[1] = j;
	  local_ga->get( vilo, vihi, &vec[0], &vld );  
          daxpyPsoci( width, value, vec, vld, hvec, vld); // we actually do not need to recompute all Hs we could transform them instead.

	}
      //vilo[0] = vlow;
      //vihi[0] = vhi; 
      vilo[1] = iroot;
      vihi[1] = iroot; 
// No need to zero out here since Hv gets reused
      local_Hga->put( vilo, vihi, &hvec[0], &vld ); // being reused as buffer space here

    }       
  /*
    At this point the vectors and H*vectors have been transformed
  */

  replaceVectors( local_Hga, num ); // pushed results back to local_ga
  return( num);
}

/* A new method that simply transforms both vec and Hvec without recomputing H*v on the final vectors.
   In this method we presume orthonormality, and simply shove the results back into
   vector, Hvector and reset state parameters
*/

//TODO double check this method for correctness
int PsociGAbasis::expandAndReplaceVectorsNew( GA::GlobalArray * g_evecs, int num )
{
#ifdef DETAILEDCHECK
  if( num > CurrentExistingRoots() ) {
    GA::error("expandAndReplaceVectorsNew: num is wrong", num);
  }
#endif

  // The plan here is to expand these eigenvectors into the existing V and Hv space, 
  // local_ga = g_array;   // the set of vectors
  // local_Hga = g_Harray
  
  int local_nbasis = CurrentNumBasis();
  int local_nsub = CurrentExistingRoots();
  
  double dZero = 0.0;
  double value;
  
  /* With this new algorithm, local_nsub should never get more than 3-4*numSoughtRoots ~ 200
     So at most we might need to fetc 100 roots * 200 subspace ~ 2,000 * 8 = 16kB of data 
     This is per core and could be alot of GA latencies. So let's simply grab them ALL and brdcst
     which is probably optimized for all machines
  */
  
  int rld, rlo[2], rhi[2];
  vector<double> root( local_nsub * num );
  // in the format root[num][local_nsub];
  
  if ( GA::nodeid() == 0 ) {
    rlo[1]=0;
    rhi[1]=local_nsub-1;
    rlo[0]=0;
    rhi[0]=num-1;
    rld = local_nsub;
    g_evecs->get( rlo,rhi, &root[0], &rld );
  }
  GA::brdcst( &root[0], local_nsub*num*sizeof(double), 0);
  
  // DO THIS IN DATA-PARALLEL
  
  int vilo[2], vihi[2];

/*
  local_ga->distribution( GA::nodeid(), vilo, vihi );
*/

  int vlow = vector_lo[0];
  int vhi  = vector_hi[0];
  int vld;
  
  if ( vector_lo[0] < 0 || vector_lo[1] < 0 ) {
     cerr << "expand NO DATA ON ME " << GA::nodeid() << endl;
     GA::Terminate();
  }
  
  int width = vhi - vlow + 1;
  vector<double> vec( width, dZero ); //move these after initial validation
  vector<double> hvec( width, dZero );

  vector<vector<double> > tvec(num, vector<double>( width,0.0) ); //move these after initial validation
  vector<vector<double> > thvec(num, vector<double>( width,0.0) ); //move these after initial validation

  int index;

  vld = 1;
  vilo[0] = vlow;
  vihi[0] = vhi;
  for(int iroot=0; iroot<num; ++iroot )
    {
      for(int j=0; j< local_nsub; ++j ) 
	{
	  index = iroot * local_nsub + j; 
	  value = root[index]; //[root][local_nsub]
	  //vilo[0] = vlow;
	  //vihi[0] = vhi;
	  vilo[1] = j; //reuse
	  vihi[1] = j;
	  local_ga->get( vilo, vihi, &vec[0], &vld );  
          daxpyPsoci( width, value, vec, vld, tvec[iroot], vld); // we actually do not need to recompute all Hs we could transform them instead.
          local_Hga->get( vilo, vihi, &hvec[0], &vld );
          daxpyPsoci( width, value, hvec, vld, thvec[iroot], vld); // we actually do not need to recompute all Hs we could transform them instead.
	}
      //vilo[0] = vlow;
      //vihi[0] = vhi; 
      vilo[1] = iroot;
      vihi[1] = iroot; 
    }       

  /*
    At this point the vectors and H*vectors have been transformed
  */

  replaceVectorsNew( tvec, thvec, num ); // pushed results back to local_ga and local_Hga

  return( num);

}

// evecs >= existing local_ga otherwise why do a compaction ?
// ONLY replace up to numSoughtRoots.
// only grab NUM vectors from the input array

int PsociGAbasis::replaceVectors( GA::GlobalArray * g_vector, int num )
{
  /*
    char * handleMessage = "PsociGAbasis:replaceVector";
    g_vector->checkHandle( handleMessage );
    
    char * handleLocalMessage = "PsociGAbasis:replaceVector";
    local_ga->checkHandle( handleLocalMessage );
  */
  
// Check for conformance..

  int ndim = scratch_ndim;

#ifdef DETAILEDCHECK
  if ( ndim > MAX_GA_VECTOR_DIMS ) {
    cerr << "ndim is " << ndim << " and cannot be greated than " << MAX_GA_VECTOR_DIMS << endl;
    GA::Terminate();
  }
#endif

/*
  int type;
  g_vector->inquire( &type, &ndim, dims );
*/
  int replace_rows = scratch_dims[0];
  int replace_cols = num; // keep (0,n-1)
  
#if 0
  if ( GA::nodeid() == 0 ) {
    cout << " current existing is " << existing_cols << endl;
    cout << " actual_roots_sought " << actual_roots_sought << endl;
  }
#endif
  
  if ( replace_rows != global_rows ) {
    cerr << "Fatal error: The nbasis for the replace data set <> that for local_ga. Was = " << replace_rows;
    cerr << " Must be " << global_rows << endl;
    GA::error("The nbasis for the replace data set <> that for local_ga", replace_rows );
  }
  
  if ( replace_cols > existing_cols || replace_cols < actual_roots_sought ) {
    cerr << "fatal error: The number of requested replace cols <  existing cols= " << replace_cols;
    cerr << " Existing size is " << existing_cols << endl;
    GA::error("The number of requested replace cols <  existing cols= ",replace_cols);
  }
  
  int basislo = 0;
  int basishi = global_rows-1; // Same for all vectors
  int lo[2], hi[2];
  //  lo[1] = -1;
  //  hi[1] = -1;
  lo[0] = basislo;
  hi[0] = basishi;
  
  int mylo[2], myhi[2];
  mylo[0] = basislo;
  myhi[0] = basishi;
  
  mylo[1] = 0;
  myhi[1] = replace_cols - 1;
  
  lo[1] = 0;
  hi[1] = replace_cols - 1;
  
  local_ga->copyPatch('N', g_vector, mylo, myhi, lo, hi );
  
  /* Since we've overwritten replace_cols from (0,replace_cols-1), we can now call 
     the RESIZE vectors methiod which wil zero out the higher order vecs, orthonormalize them and compute
     a fresh dset of H*v.
  */ 
  int numr = resizeVectors( replace_cols );
  return( numr );
}

/* Here we simply shove the result back into GA space startting at 0
*/
int PsociGAbasis::replaceVectorsNew( vector<vector<double> >&tvec, vector<vector<double> >&thvec, int num )
{

  /*
    char * handleMessage = "PsociGAbasis:replaceVector";
    g_vector->checkHandle( handleMessage );
    
    char * handleLocalMessage = "PsociGAbasis:replaceVector";
    local_ga->checkHandle( handleLocalMessage );
  */
  
// Check for conformance..

  int vlow = vector_lo[0]; // local parts of vec and Hvec
  int vhi  = vector_hi[0];

  if ( tvec.size() != thvec.size() || tvec.size() != num ) GA::error(" replaceVectorsNew failed",tvec.size() );

#if 0
  if ( GA::nodeid() == 0 ) {
    cout << "new current existing is " << existing_cols << endl;
    cout << "new actual_roots_sought " << actual_roots_sought << endl;
  }
#endif
  
  int mylo[2], myhi[2];
  mylo[0] = vlow;
  myhi[0] = vhi;
  
  int n = 1;

/* We only really need to do this for the Scatter based techniques 
   Still much faster to transform the products

   After this we upload as data parallel
*/
  local_Hga->zero();

  for(int iroot=0; iroot<num; ++iroot )
  {

  mylo[1] = iroot;
  myhi[1] = iroot; 
  
  local_Hga->put( mylo, myhi, &thvec[iroot][0], &n );
  local_ga->put( mylo, myhi, &tvec[iroot][0], &n );
  }
  
/* skip orthonormalization and don't worry about zeroing out higher valued GA columns
*/

 present_matrixProducts = num; // start with transformed products
 required_matrixProducts = 0; // Restart subspace >= soughtNum
 existing_cols = num;

 return( num );
}


// We do NOT update actual_roots_sought since this is intended for use by subspace iteration which grows/shrinks in time.
// Appends vectors stored in g_vector to the vector space
// We still need to update H*v, though

int PsociGAbasis::appendVectors( GA::GlobalArray * g_vector )
{

#ifdef DETAILECHECK
    char * handleMessage = "PsociGAbasis:appendVector";
    g_vector->checkHandle( handleMessage );
    
    char * handleLocalMessage = "PsociGAbasis:appendVector";
    local_ga->checkHandle( handleLocalMessage );
#endif

  // Check for conformance..
  
    int ndim = g_vector->ndim();

#ifdef DETAILEDCHECK
    if ( ndim > MAX_GA_VECTOR_DIMS ) {
      cerr << "ndim is " << ndim << " and cannot be greated than " << MAX_GA_VECTOR_DIMS << endl;
      GA::Terminate();
    }
#endif

/*
  int type;
  g_vector->inquire( &type, &ndim, dims );
*/

  int append_cols = scratch_dims[1];

#if 0
  if ( GA::nodeid() == 0 ) {
    cout << " current existing is " << existing_cols << endl;
    cout << " desired cols is " << desired_cols << endl;
    cout << " actual_roots_sought " << actual_roots_sought << endl;
    cout << " append_cols " << append_cols << endl;
  }
#endif

#ifdef DETAILEDCHECK
//  if ( append_rows != global_rows ) {
//    cerr << "Fatal error: The nbasis for the append data set <> that for local_ga. Was = " << append_rows;
//    cerr << " Must be " << global_rows << endl;
//    GA::Terminate();
//  }
//  desired_cols = existing_cols + append_cols;
//  if ( desired_cols > global_cols ) {
//    cerr << "fatal error: The number of requested append cols exceeds global data size Would be = " << desired_cols;
//    cerr << " Maximum size is " << global_cols << endl;
//    GA::Terminate();
//  }
#endif
  
  int basislo = 0;
  int basishi = global_rows-1; // Same for all vectors
  int lo[2], hi[2];
  lo[1] = -1;
  hi[1] = -1;
  lo[0] = basislo;
  hi[0] = basishi;
  
  int mylo[2], myhi[2];
  mylo[0] = basislo;
  myhi[0] = basishi;
  
  int global_index;
  for (int i = 0; i < append_cols; ++i ) {
    global_index = existing_cols + i;
    lo[1] = global_index;
    hi[1] = global_index;
   
    mylo[1] = i;
    myhi[1] = i;
    
#if 0
    cout << " LO " << lo[0] << " " << hi[0] << lo[1] << " " << hi[1] << endl;
    cout << " MYLO " << mylo[0] << " " << myhi[0] << mylo[1] << " " << myhi[1] << endl;
#endif
    local_ga->copyPatch('N', g_vector, mylo, myhi, lo, hi );
  }
  /* older version for ortho ALL vectors 
     int or_status = orthonormalizeAllVectors();
  */
  
#ifdef NEWORTHO
  int or_status = orthonormalizeAppendedVectors(existing_cols, append_cols);
#else
  int or_status = orthonormalizeAppendedVectorsGA(existing_cols, append_cols);
#endif

  existing_cols += append_cols; //update index
  
  /* Maybe we should not do this automatically
     ONLY update appended vectors worth of H*v
  */
  
  required_matrixProducts = existing_cols - present_matrixProducts;

#if 0
   vector<double> rvec( global_rows, 0.0 );
   int ld=1;
   lo[0]=0;
   lo[1]=1;
   hi[0]=global_rows-1;
   hi[1]=1;
   local_ga->get( lo, hi, &rvec[0], &ld);
   cout << "DUMP RVEC AFTER NORM COEFS for 1" << endl;
   for(int k=0; k< global_rows; ++k ){
   cout << setprecision(15) << rvec[k] << " " ;
   }
   cout << endl;
#endif
  
   int numH = computeHvProducts( present_matrixProducts, required_matrixProducts );
   present_matrixProducts = existing_cols;
   required_matrixProducts = 0;
   
   return( or_status );
}

//Collective
void PsociGAbasis::getAllNorms( vector<double> & norm, pair<int,double> & time )
{
  int rank = GA::nodeid();
  double timein = psociTime();
  if ( PsociGAbasis::getAllNorms( norm ) != 0 ) {
    cerr << "getAllNorms-2 failed with error: Aborting " << endl;
    GA::Terminate();
  }
  time.first = rank;
  time.second = psociTime() - timein;
}

/* We could copy GA  matrices into a set of GA_vectors and simply perform the norm1. But I am
   betting the time will be significantly influenced by the GA creates that would be needed. So
   Do the copy patch approach
   
   On output returns the set of Frobenius Norms[existing_rows];
*/
// Collective: 
int PsociGAbasis::getAllNorms( vector<double> & norms )
{
  int normstatus = 0;
  int rank = GA::nodeid();
  const double precision = 0.01 * SQ_MIN_PRECISION;
  
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:getAllNorms";
  local_ga->checkHandle( handleMessage );
  
    if ( norms.size() > 0 ) {
    cout << "getAllNorms: Warning: expected an empty norms array. Emptying and proceeding " << norms.size() << endl;
    }
#endif

    norms.resize( existing_cols, 0.0 );
    
  // Will place temporary values in the local array double V[roots] then global sum the whole thing
  // Or we could stick results into a GA vector of length nroots ?
  // Only process existing_rows. If/when new vectors are added to the groups, then update existing to reflect this.
    
  // Note distribution provides dims in C-indexing style (0,n-1)
    
    int lo[2], hi[2];
    int rootlo = 0;
    int roothi = existing_cols-1;
    int basislo = 0;
    int basishi = global_rows-1;
    
    double value;
#if 0
    cout << "rank is " << rank << " lo 0 and hi 0 are " << rootlo << " " << roothi << endl;
    cout << "rank is " << rank << " lo 1 and hi 1 are " << basislo << " " << basishi << endl;
#endif
    
    lo[0] = basislo;
    hi[0] = basishi;
    
    for (int i=rootlo; i <= roothi; ++i ) 
      {
	lo[1]=i; // Process i-root at a time
	hi[1]=i;      
	value = sqrt(local_ga->ddotPatch('N', lo, hi, local_ga, 'N', lo, hi ));
	norms[i] = value;
	if ( value < precision ) {
	  cout << "Warning: getAllnodrms: An inordinately low Frobenius norm " << value;
          cout << " For root " << i+1<<" Fortran-style "  << " precision is " << precision << endl;
	  normstatus = 1;
	}
      }
    return( normstatus );
}

//Collective: 
/* 
   start is an absolute number of roots
*/
//TODO Need to check the use of precision terms. They are getting abit out of sync
int PsociGAbasis::getAppendedNorms( vector<double> & norms, int start, int num )
{
  int normstatus = 0;
  const double precision = 0.01 * SQ_MIN_PRECISION;
  
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:getAppendedNorms";
  local_ga->checkHandle( handleMessage );
  if ( norms.size() > 0 ) {
    cout << "getAllNorms: Warning: expected an empty norms array. Emptying and proceeding " << norms.size() << endl;
  }
#endif
  
#ifdef DETAILEDCHECK
  norms.clear();
#endif
  norms.resize( num, 0.0 );
  
  int lo[2], hi[2];
  int rootlo = start;
  int basislo = 0;
  int basishi = global_rows-1;
  
  double value;
#if 0
  int rank = GA::nodeid();
  cout << "rank is " << rank << " lo 0 and hi 0 are " << rootlo << " " << roothi << endl;
  cout << "rank is " << rank << " lo 1 and hi 1 are " << basislo << " " << basishi << endl;
#endif
  
  lo[0] = basislo;
  hi[0] = basishi;
  int i;
  for (int ind=0; ind < num; ++ind )
    {
      i = ind + rootlo;
      lo[1]=i; // Process i-root at a time
      hi[1]=i;      
      value = sqrt(local_ga->ddotPatch('N', lo, hi, local_ga, 'N', lo, hi ));
      norms[ind] = value;
      if ( value < precision ) {
	cout << "getAppended norms: Warning: An inordinately low Frobenius norm " << value << endl;
        cout << " For root " << i+1<<" Fortran-style "  << " precision is " << precision << endl;

	normstatus = 1;
      }
    }
  return( normstatus );
}
//TODO Need to check the use of precision terms. They are getting abit out of sync
int PsociGAbasis::getAppendedNormsDataParallel( vector<double> & norms, int start, int num )
{
  int normstatus = 0;
// Currently SQ_MIN_PRECISION is set to zero.
  
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:getAppendedNormsDataParallel";
  local_ga->checkHandle( handleMessage );
  if ( norms.size() > 0 ) {
    cout << "getAllNormsDataParallel: Warning: expected an empty norms array. Emptying and proceeding " << norms.size() << endl;
  }
#endif
  
  const double dZero = 0.0;
  norms.clear();
  norms.resize( num, 0.0 );
  
  int lo[2], hi[2];
  int rootlo = start;
  
  lo[0] = vector_lo[0];
  hi[0] = vector_hi[0];
  int width = hi[0] - lo[0] + 1;
  vector<double> vec( width, dZero );

  if ( lo[0] < 0 ) cout << "BAD DISTRIBUTION" << endl;

  double value;
#if 0
  cout << "rank is " << rank << " lo 0 and hi 0 are " << rootlo << " " << roothi << endl;
  cout << "rank is " << rank << " lo 1 and hi 1 are " << basislo << " " << basishi << endl;
#endif
  
  int i;
  int ld = 1;
  char op[] = "+";
  int numElems  = 1;
  
  for (int ind=0; ind < num; ++ind )
    {
      i = ind + rootlo;
      lo[1]=i; // Process i-root at a time
      hi[1]=i;      
      local_ga->get( lo, hi, &vec[0], &ld );

      value = ddotPsoci( width, vec, ld, vec, ld );
      GA::gop( (double *) &value, numElems, op );

      norms[ind] = sqrt( value );
//      if ( value < precision ) {
//	cout << "Appended norms: data parallel: Warning: An inordinately low Frobenius norm " << value << endl;
//        cout << " For root " << i+1<<" Fortran-style "  << "value is " << value<< " precision is " << precision << endl;
//	normstatus = 1;
//      }
    }
  return( normstatus );
}

// Collective - timer wrapper
void PsociGAbasis::normAllVectors( pair<int,double> & time )
{
  int rank = GA::nodeid();
  double timein = psociTime();  
  PsociGAbasis::normAllVectors( );
  time.first = rank;
  time.second = psociTime() - timein;
}

//Collective - normalizes the GA in-place
//Only process against the set of existing rows......
//Presumes a unit metric and does not check for orthogonality ( use something else for that )

void PsociGAbasis::normAllVectors( )
{
  const double precision = 0.01 * SQ_MIN_PRECISION;

#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:normAllVecs";
  local_ga->checkHandle( handleMessage );
#endif
  
  vector<double> norms;
  int normstatus = getAllNorms( norms );
  if ( normstatus != 0 ) {
    cerr << " A vector was to within " << precision << " of zero. They should all be checked/normalized/projected " << endl;
  }
  int lo[2], hi[2];
  
  //Basis size
  lo[0]=0;
  hi[0]=global_cols-1;
  
  // Process all existing_rows
  for (int i=0; i < existing_cols; ++i ) 
    {
      lo[1]=i;
      hi[1]=i;
      if ( norms[i] >= precision ) {
        double scale = 1.0/norms[i];
	local_ga->scalePatch( lo, hi, &scale );
      } else {
	cerr << "Error: norms[i] close to zero: You may need to orthogonalize the vectors and project out dependants " << endl; 
        cerr << "Current erroneous norm is " << norms[i] << " at root number (i+1) " << i+1 << endl;
	GA::Terminate();
      }
    }
}

//Remove GA-based ddot in favor of more explicit control over memory use. 
//Collective are any vectors NOT orthogonal. Compute overlap matrix on I^basis and report number of non-zeros

int PsociGAbasis::checkAllOverlap( vector<pair<int,int> > & overlap )
{
  const double precision = 0.01 * SQ_MIN_PRECISION;
  int orthstatus = 0;

#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:checkAllOverlap";
  local_ga->checkHandle( handleMessage );
#endif
  
  // both are maxsef X roots. distributon of maxsef must be the same
/*
  int vilo[2], vihi[2];
  local_ga->distribution( GA::nodeid(), vilo, vihi );
*/
  int vlow = vector_lo[0];
  int vhi  = vector_hi[0];
  int width = vhi - vlow + 1;
  
  int roothi = existing_cols;
  
  int rank = GA::nodeid();
  double value;
  double dZero = 0.0;
  
  int lo[2], hi[2];
  
  lo[0] = vlow;
  hi[0] = vhi;
  int ld = 1;
  
  vector<double> veci( width, dZero );
  vector<double> vecj( width, dZero );
  
  for (int i = 0; i < roothi; ++i )
    { 
      lo[1] = i;
      hi[1] = i;
      local_ga->get( lo, hi, &veci[0], &ld );
      
      for( int j = 0; j < i; ++j )
        { 
	  lo[1] = j;
	  hi[1] = j;
	  local_ga->get( lo, hi, &vecj[0], &ld );
          value = ddotPsoci( width, veci, ld, vecj, ld );
#if 0
	  cout << "Overlap local ddot rank " << rank << " I+1 is " << i+1 <<  " J+1 is " << j+1 << " value " << value << endl;
#endif
	  if ( value >= precision ) {
	    overlap.push_back( pair<int,int>(i,j) ); 
	    // cout << "Warning: I and J are non-zero to greater than precision" << value << endl;
	    orthstatus = 1;
	  }
	}
    }
  return( orthstatus );
}

//Collective are any vectors NOT orthogonal. Compute overlap matrix on I^basis and report number of non-zeros
int PsociGAbasis::checkAllOverlapGA( vector<pair<int,int> > & overlap )
{

  const double precision = 0.01 * SQ_MIN_PRECISION;
  int orthstatus = 0;

#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:checkAllOverlap";
  local_ga->checkHandle( handleMessage );
#endif
  
  int lo[2], hi[2]; // Vector number I 
  int lo2[2], hi2[2]; //Vector number J
  
  int roothi = existing_cols;
  
  int rank = GA::nodeid();
  double value;
  
  lo[0] = 0;
  hi[0] = global_rows-1;
  lo2[0] = 0;
  hi2[0] = global_rows-1; //Full vector lenghts for I and J
  
  for (int i = 0; i < roothi; ++i )
    { 
      lo[1]=i; // Process i-root 
      hi[1]=i;      
      for( int j = 0; j < i; ++j )
        { 
	  lo2[1] = j;
	  hi2[1] = j; // Process j-root;
	  value = sqrt(local_ga->ddotPatch('N', lo, hi, local_ga, 'N', lo2, hi2 ));
#if 0
	  cout << "Overlap rank " << rank << " I+1 is " << i+1 <<  " J+1 is " << j+1 << " value " << value << endl;
#endif
	  if ( value >= precision ) {
	    overlap.push_back( pair<int,int>(i,j) ); 
	    // cout << "Warning: I and J are non-zero to greater than precision" << value << endl;
	    orthstatus = 1;
	  }
	}
    }
  return( orthstatus );
}

int PsociGAbasis::orthonormalizeAllVectors( pair<int,double> & time )
{
  int rank = GA::nodeid();
  double timein = psociTime();
  int status = orthonormalizeAllVectors();
  time.first = rank;
  time.second = psociTime() - timein;
  return( status );
}

//Now we need a modified gram-schmidt orthogonalization 
/** This routine takes on input the set of existing_rows and mutually orthogonalizes them all
    This is probably not necessary all the time HOWEVER, it does result in good precision. 
    ALternative orthogonalizers processing only augmented vectors are also available......
    
    This will trigger a rebuild of all H*v products in the new subspace.
*/

// Do the normalization twice to retain precision
// return 0 on success or k for failure of the kth vector

int PsociGAbasis::orthonormalizeAllVectors()
{
  int orthostatus = 0;
  int rank = GA::nodeid();
  const double precision = 0.01 * SQ_MIN_PRECISION;
  
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:orthonormalizeAllVectors";
  local_ga->checkHandle( handleMessage );
#endif
  
  int lo[2], hi[2]; // I vector
  int lo2[2], hi2[2]; // J vector
  
  //Set up basic Patch coordinates for the daxpy
  int basishi = global_rows - 1;
  int rowshi = existing_cols - 1;

#if 0
  cout << "ortho: rank is " << rank << "existing roots is " << existing_rows << endl;
#endif
  
  // Perform path-GA modified gram-schmidt orthogonalization (QR) step followed by normalization step
  // Compute new orthogonal basis
  
  lo[0] = 0;
  hi[0] = basishi; //Al of I
  lo2[0] = 0;
  hi2[0] = basishi; //All of J 
  
  double proj;
  double norm;
  double overlap;
  double dONE = 1.0;
  double dNegate = -1.0;
  
  for(int iter=0; iter < NUM_ORTHO_REPS ; ++iter ) // Do it twice
    {
      for (int i = 0; i <= rowshi; ++i ) {
	lo[1] = i;
	hi[1] = i;
	for (int j = 0; j < i ; ++j ) {
	  lo2[1] = j;
	  hi2[1] = j;
	  /* Projection operator:  proj_i(j)
	     Vector(i) = Vector(i) + (-1) * overlap * Vector(j) . No need to normalize overlap as that is handled 
	     implicitly
	  */
	  overlap = (dNegate) * local_ga->ddotPatch('N', lo, hi, local_ga, 'N', lo2, hi2 );
	  norm = sqrt(local_ga->ddotPatch('N', lo2, hi2, local_ga, 'N', lo2, hi2 ));
	  
	  if ( norm > precision ) {
	    proj = overlap / norm; 
	    //     if ( abs(overlap) >= sq_precision ) local_ga->addPatch(&dONE, local_ga, lo, hi, &proj, local_ga, lo2, hi2, lo, hi );
	    local_ga->addPatch(&dONE, local_ga, lo, hi, &proj, local_ga, lo2, hi2, lo, hi );
	  } else {
	    cerr << "Error orthogonalizingAll Probably  a linear dependency: vector i = " << i+1 << " Should be removed " << endl;
	    GA::error("Error orthogonalizingAll:linear dep ", i+1);
	  }   
#if 0
	  cout << "OVERLAP (-1) between i+1 and j+1 = " << i+1 << " " << j+1 << " is " << overlap << endl;
	  cout << "Partial updated matrix print " << endl;
	  local_ga->print();
#endif
	  
	}
      }  
    }
  return( orthostatus );
}

/** In this case we want to save time by simply orthogonalizing only
    the new vectors. This is a cheap operation BUT we did not want to 
    orthoginalize all vectors and either recompute all H*v and allow them
    to become rotated relative to one another ( even thoought that is unitary) 
    
    Othonormaliuze starting at "start" and ending at "< num+start")
*/

// Do the normalization twice to retain precision
// return 0 on success or k for failure of the kth vector
// Validated

int PsociGAbasis::orthonormalizeAppendedVectorsGA(int start, int num )
{
  int orthostatus = 0;
  int rank = GA::nodeid();
  const double precision = MIN_PRECISION;
  
  
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:orthonormalizeAppendedVectorsGA";
  local_ga->checkHandle( handleMessage );
#endif
  
  int lo[2], hi[2]; // I vector
  int lo2[2], hi2[2]; // J vector
  
  int basishi = global_rows;
  
#if 0
  cout << "ortho:set:: rank is " << rank << "existing roots is " << existing_rows << endl;
#endif
  
  // Perform path-GA modified gram-schmidt orthogonalization (QR) step followed by normalization step
  // Compute new orthogonal basis
  
  lo[0] = 0;
  hi[0] = basishi-1; //Al of I
  lo2[0] = 0;
  hi2[0] = basishi-1; //All of J 
  
  double norm;
  double overlap;
  double proj;
  double dONE = 1.0;
  double dNegate = -1.0;
  
/* orthogonalize i against the set j
*/
  for(int iter=0; iter<NUM_ORTHO_REPS; ++iter)  // Do it twice for precision
    {
      for (int i = start; i < start+num; ++i ) {
	lo[1] = i;
	hi[1] = i;
	for (int j = 0; j < i ; ++j ) {
	  lo2[1] = j;
	  hi2[1] = j;
	  /* Projection operator:  proj_i(j)
	     Vector(i) = Vector(i) + (-1) * overlap * Vector(j) . No need to normalize overlap as that is handled 
	     implicitly
	  */
	  overlap = (dNegate) * local_ga->ddotPatch('N', lo, hi, local_ga, 'N', lo2, hi2 );
	  norm = sqrt(local_ga->ddotPatch('N', lo2, hi2, local_ga, 'N', lo2, hi2 ));
	  
	  if ( norm > precision ) {
	    proj = overlap / norm;
	    //     if ( abs(overlap) >= sq_precision ) local_ga->addPatch(&dONE, local_ga, lo, hi, &proj, local_ga, lo2, hi2, lo, hi );
	    local_ga->addPatch(&dONE, local_ga, lo, hi, &proj, local_ga, lo2, hi2, lo, hi );
	  } else {
	    cerr << "Error orthogonalizingAppend: Probably  a linear dependency: vector i = " << i+1 << " Should be removed " << endl;
	    GA::error("Error orthogonalizingSet:linear dep ", i+1);
	  } 
	}     
      }  
    }
  normAppendedVectors( start, num );
  return( orthostatus );
}

/** In this case we want to save time by simply orthogonalizing only
    the new vectors. This is a cheap operation BUT we did not want to 
    orthoginalize all vectors and either recompute all H*v and allow them
    to become rotated relative to one another ( even thoought that is unitary) 

    Othonormalize starting at "start" and ending at "< num+start")
*/

// Do the normalization twice to retain precision
// return 0 on success or k for failure of the kth vector

int PsociGAbasis::orthonormalizeAppendedVectors(int start, int num )
{
  int orthostatus = 0;
//  int rank = GA::nodeid();
  const double precision = MIN_PRECISION;
  
#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:orthonormalizeAppendedVectors";
  local_ga->checkHandle( handleMessage );
#endif
  
  int lo[2], hi[2]; // I vector
  int lo2[2], hi2[2]; // J vector

  //  int basishi = global_rows;
  //  int rowshi = existing_cols;

#if 0
  cout << "ortho:set:: rank is " << rank << "existing roots is " << existing_rows << endl;
#endif
  
  // Perform path-GA modified gram-schmidt orthogonalization (QR) step followed by normalization step
  // Compute new orthogonal basis

//  int vilo[2], vihi[2];
/*
  local_ga->distribution( GA::nodeid(), vilo, vihi );
*/
  int vlow = vector_lo[0];
  int vhi  = vector_hi[0];
  int ld=1;
  
  lo[0] = vlow;
  hi[0] = vhi; //Al of I
  lo2[0] = vlow;
  hi2[0] = vhi; //All of J 
  
  int width = vhi - vlow + 1;
  
  double norm;
  double overlap;
  double proj;
  double dZero = 0.0;
  double dNegate = -1.0;

  double pair[2]; //local overlap/norm
 
  char op[] = "+";
  int numElems  = 2;
  
  /* orthogonalize i against the set j
   */

  vector<double> veci( width, dZero );
  vector<double> vecj( width, dZero );

  for(int iter=0; iter<NUM_ORTHO_REPS; ++iter)  // Do it NUM_ORTHO_REPS for precision
    {
//      vector<double> veci( width, dZero );
      vector<double> vecj( width, dZero );
      
      for (int i = start; i < start+num; ++i ) {
	lo[1] = i;
	hi[1] = i;
	local_ga->get( lo, hi, &veci[0], &ld );
	
	for (int j = 0; j < i ; ++j ) {
	  lo2[1] = j;
	  hi2[1] = j;
	  local_ga->get( lo2, hi2, &vecj[0], &ld );
	  
	  /* Projection operator:  proj_i(j)
	     Vector(i) = Vector(i) + (-1) * overlap * Vector(j) . No need to normalize overlap as that is handled 
	     implicitly
	  */
	  
          pair[0] = (dNegate) * ddotPsoci( width, veci, ld, vecj, ld );
          pair[1] = ddotPsoci( width, vecj, ld, vecj, ld );
	  GA::gop( (double *) &pair, numElems, op );

	  overlap =  pair[0] ;
	  norm = sqrt ( pair[1] );

	  if ( norm > precision ) {
	    proj = overlap / norm;
	    //    if ( abs(overlap) >= sq_precision ) local_ga->addPatch(&dONE, local_ga, lo, hi, &proj, local_ga, lo2, hi2, lo, hi );
	    //    local_ga->addPatch(&dONE, local_ga, lo, hi, &proj, local_ga, lo2, hi2, lo, hi );    
          daxpyPsoci( width, proj, vecj, ld, veci, ld );

	  } else {
	    cerr << "Error new orthogonalizingAppend: Probably  a linear dependency: vector i = " << i+1 << " Should be removed " << endl;
	    GA::error("Error new orthogonalizingSet:linear dep ", i+1);
	  } 
	}     
	local_ga->put(lo, hi, &veci[0], &ld );
      }  
    }  
  normAppendedVectors( start, num );
  return( orthostatus );
}


/* Normalize num vectors starting at start.
   For now go ahead and get all norms. but we do need to prune down that list eventually
*/
void PsociGAbasis::normAppendedVectors(int start, int num )
{
  //int rank = GA::nodeid();
  //const double precision = MIN_PRECISION;
  const double precision = 0.0;
  
#ifdef DETAILEDCHECK
    char * handleMessage = "PsociGAbasis:normAppendedVectors";
    local_ga->checkHandle( handleMessage );
#endif
  
  int lo[2], hi[2]; // I vector
  
  vector<double> norms;
  
  //Basis size
  lo[0]=0;
  hi[0]=global_rows-1;
  
  double dOne= 1.0;
  
  norms.clear();
#ifdef NEWGETAPPEND
  int normstatus = getAppendedNormsDataParallel(norms,start,num );
#else
  int normstatus = getAppendedNorms(norms,start,num );
#endif
  
  if ( normstatus != 0 ) {
    cout << "normAppendedVectors: A vector was to within " << precision << " of zero. They should all be checked/normalized/projected " << endl;
    cout << " Current range is  " << start<< " "<< num << endl; 
  }
  
#if 0
    for(int i=0; i< num; ++i ) {
    cout << "NORMS fetched are " << norms[i] << endl;
    }
#endif
  
    // Process all existing_rows
    int i;
    for (int ind=0; ind < num; ++ind )
      {
	i = ind + start;
	lo[1]=i;
	hi[1]=i;
	if ( norms[ind] > precision ) {
	  double scale = dOne/norms[ind];
	  local_ga->scalePatch( lo, hi, &scale );
	} else {
	  cerr << "normAppendedVectors: Error: norms[i] close to zero: You may need to orthogonalize the vectors and project out dependants " << endl;
	  cerr << "normAppendedVectors: Current erroneous norm is " << norms[ind] << " at root number (i+1) " << i+1 << " " << norms[ind] << endl;
	  GA::error( "normAppendedVectors: Current erroneous norm is ", i+1);
	}
      }
}
/** This method is called when attempting to append a vector fails because of 
    exceeding maxsub.  In this case, we want to simply remove all non-converged vectors
    (zero), and reset existing_cols. Only keep converged_cols number of vectors untouched
*/

void PsociGAbasis::compactVectors()
{
  int highest_converged_root = converged_cols;
  int highest_root = existing_cols;
  
  // Set up patch dims

#ifdef DETAILEDCHECK
  char * handleMessage = "PsociGAbasis:compact";
  local_ga->checkHandle( handleMessage );
  
  char * handleHvMessage = "PsociGAbasis:compact:Hv";
  local_Hga->checkHandle( handleHvMessage );
#endif
  
  int dims[MAX_GA_VECTOR_DIMS];
  int ndim = local_ga->ndim();

#ifdef DETAILEDCHECK
  if ( ndim > MAX_GA_VECTOR_DIMS ) {
    cerr << "ndim is " << ndim << " and cannot be greated than " << MAX_GA_VECTOR_DIMS << endl;
    GA::Terminate();
  }
#endif

//  int type;
//  local_ga->inquire( &type, &ndim, dims );

  //global_rows = dims[0];
  //global_cols = dims[1];
  if ( highest_root > vector_dims[1] ) {
    cerr << "Inconsistent array dimensions " << highest_root << " " << vector_dims[1] << endl;
    GA::Terminate();
  }
  
  if ( highest_converged_root == 0 ) {
    cerr << "highest_converged_root == 0. This is unexpected " << endl;
    GA::Terminate();
  }
  
  // TODO MAYBE WE CAN DO THIS
  //  local_ga->zeroPatch( lo, hi ); 
  //  local_Hga->zeroPatch( lo, hi );
  
  existing_cols = highest_converged_root;
  desired_cols = existing_cols;
}


/** Not a real random number generator BUT, it is only being used to 
    generate guesses for augmented CI vectors which will be subsequently
    orthonormalized.
    On return: a double between 0 and 10
*/
void PsociGAbasis::fetchPseudoRandom( int number, int seed, vector<double> & buf )
{
  double dTen = 10.0;
  double dLow = 0.0;
  srand( seed );
  double random_double; // Between 0 and 10
  for (int i = 0; i < number; ++i ) {
    random_double =  rand() * dTen / RAND_MAX + dLow; 
    //cout << ++list << " random double " << random_double << endl;
    buf.push_back( random_double );
  }
}


/** Compute the H*v products for the vectors(start:end) beginning at start and ending at end-1.
    local_ga, and local_Hga

    Need to account for precision at some point but for now use full precision
*/
int PsociGAbasis::computeHvProducts( int start, int num )
{
  int end = start + num;
//  int dims[2];
//  int ndim = 2;
//  int type;
//  local_ga->inquire( &type, &ndim, dims );
  
#ifdef DETAILEDCHECK
  if ( vector_dims[1] < num ) {
    cerr << "Failed computeHvProducts: num too large :" << num << endl;
    GA::error("Failed computeHvProducts: num too large :", vector_dims[1] );
  }
#endif
  //Assume g_Hv distribution is okay
  //Need to pass in dims/lo/hi for perforance reasons

// SELECT the style of MatrixVector product to use
// CHOICES -DGATHERHV ( bounded) -DMATRIXPRODUCTGOP (fast - smaller problems), Old style default exists but never used
// CHOICE options also exists for -DSUPPLEMENTAL (split GA space) 

/*
#ifdef GATHERHV
#ifdef SUPPLEMENTAL 
 int vcount = local_hamiltonian->matrixVectorProductsIncoreGOPSplitGatherScatter( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
#else
//int vcount = local_hamiltonian->matrixVectorProductsIncoreGOPSplit( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
#endif
#endif

#ifdef MATRIXPRODUCTGOP
#ifdef GATHERHV
    GA::error("Both GATHERHV and MATRIXPRODUCTGOP specified. CHoose one",1);
#endif
#ifndef SUPPLEMENTAL 
  GA::error("To use SUPPLEMENTAL you must also preprocess with -DMATRIXPRODUCTGOP",-1);
#endif
//int vcount = local_hamiltonian->matrixVectorProductsIncoreGOPSplit( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
#endif
*/

//if ( GA::nodeid() == 0 ) cout << "gather " << endl;
//int vcount = local_hamiltonian->matrixVectorProductsIncoreGOPSplitGatherScatter( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );


//if ( GA::nodeid() == 0 ) cout << " GOP version " << endl;

#ifdef SQUAREHAMILTONIAN
GA::error(" NO reason to use SQUAREHAMILTONIAN any more",1);
int vcount = local_hamiltonian->matrixVectorProductsIncoreGOPSplitChunk(start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );


#else

//if ( GA::nodeid() == 0 ) cout << " normal GOPSplit" << endl;
//int vcount = local_hamiltonian->matrixVectorProductsIncoreGOPSplit( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

/* THis is pretty good
if ( GA::nodeid() == 0 ) cout << " NEW Chunk fetch replicated accume " << start <<" "<<endl<<endl;

*/
/* best
*/

int vcount;
if ( l_scatter_products ) {

#ifndef DIRECTMATRIXPRODUCTS
if ( GA::nodeid() == 0 ) cout << " matrixVectorProductsSplitChunkScatterAccumulate new ACCUMULATE " << start <<" "<<endl<<endl;
 vcount = local_hamiltonian->matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
#endif
#ifdef DIRECTMATRIXPRODUCTS
//if ( GA::nodeid() == 0 ) cout << "DIRECT experimental  matrixVectorProductsSplitChunkScatterAccumulateDirect new ACCUMULATE " << start <<" "<<endl<<endl;
//int vcount = local_hamiltonian->matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirect( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

if ( GA::nodeid() == 0 ) cout << "DIRECT + FETCH experimental  matrixVectorProductsSplitChunkScatterAccumulateDirectFetch new ACCUMULATE " << start <<" "<<endl<<endl;
 vcount = local_hamiltonian->matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirectFetch( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
#endif

#ifdef INNERSCATTER
if ( GA::nodeid() == 0 ) cout << " using INNERSCATTER " << endl;
#endif

} else {

/* a this point we could change the value of l_nxcheck to fit memory constraints */

 vcount = local_hamiltonian->matrixVectorProductsSplitChunkMultipleJchunks( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

//if ( GA::nodeid() == 0 ) cout << " vcount is " << vcount<<" "<<start<<" "<<end<<endl;

}


#ifdef  NEWACCUMULATEDSCATTER
//if ( GA::nodeid() == 0 ) cout << "NEW SCATTER  Total memory per core for M*v is (MB)" << local_hamiltonian->matrixVectorProductsMemoryPerCore()/1000000.0 << endl;
//if ( GA::nodeid() == 0 ) cout << " matrixVectorProductsSplitChunkScatterAccumulate" << start <<" "<<endl<<endl;
//int vcount = local_hamiltonian->matrixVectorProductsSplitChunkScatterAccumulate( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

// New version calling lower level GA scatter
//if ( GA::nodeid() == 0 ) cout << "NEW LOWER LEVEL SCATTER  Total memory per core for M*v is (MB)" << local_hamiltonian->matrixVectorProductsMemoryPerCore()/1000000.0 << endl;
//if ( GA::nodeid() == 0 ) cout << " matrixVectorProductsSplitChunkScatterAccumulate " << start <<" "<<endl<<endl;
//int vcount = local_hamiltonian->matrixVectorProductsSplitChunkScatterAccumulateLowLevel( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

#ifndef DIRECTMATRIXPRODUCTS
//keepif ( GA::nodeid() == 0 ) cout << " matrixVectorProductsSplitChunkScatterAccumulate new ACCUMULATE " << start <<" "<<endl<<endl;
//keepint vcount = local_hamiltonian->matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
#endif

#ifdef DIRECTMATRIXPRODUCTS
//if ( GA::nodeid() == 0 ) cout << "DIRECT experimental  matrixVectorProductsSplitChunkScatterAccumulateDirect new ACCUMULATE " << start <<" "<<endl<<endl;
//int vcount = local_hamiltonian->matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirect( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

//keepif ( GA::nodeid() == 0 ) cout << "DIRECT + FETCH experimental  matrixVectorProductsSplitChunkScatterAccumulateDirectFetch new ACCUMULATE " << start <<" "<<endl<<endl;
//keepint vcount = local_hamiltonian->matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirectFetch( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
#endif

#ifdef INNERSCATTER
//keepif ( GA::nodeid() == 0 ) cout << " using INNERSCATTER " << endl;
#endif


#else
//keepif ( GA::nodeid() == 0 ) cout << " matrixVectorProductsSplitChunkRepAccumulate " << start <<" "<<endl<<endl;
//keepif ( GA::nodeid() == 0 ) cout << " Total memory per core for M*v is (MB)" << local_hamiltonian->matrixVectorProductsMemoryPerCore()/1000000.0 << endl;
//keepint vcount = local_hamiltonian->matrixVectorProductsSplitChunkRepAccumulate( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
#endif

// DO NOT USE BELOW HERE
/* OLD
if ( GA::nodeid() == 0 ) cout << " NEW LOCAL VECTOR fetch replicated accume " << start <<" "<<endl<<endl;
if ( GA::nodeid() == 0 ) cout << " Total memory per core for M*v is (MB)" << local_hamiltonian->matrixVectorProductsMemoryPerCore()/1000000.0 << endl;
int vcount = local_hamiltonian->matrixVectorProductsSplitLocalVectorRepAccumulate( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );
*/

// KEEP THIS
//if ( GA::nodeid() == 0 ) cout << " gatherscatter " << start <<" "<<endl<<endl;
//int vcount = local_hamiltonian->matrixVectorProductsIncoreGatherScatter( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

// expeirmental
//if ( GA::nodeid() == 0 ) cout << "block vector GATHER  gatherscatter " << start <<" "<<endl<<endl;
//int vcount = local_hamiltonian->matrixVectorProductsIncoreGatherScatterNoVectorGet( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

//if ( GA::nodeid() == 0 ) cout << " old gatherscvatter " << endl;
//int vcount = local_hamiltonian->matrixVectorProductsIncoreGOPSplitGatherScatter( start, end, local_ga, &vector_dims[0], local_Hga, &vector_dims[0],&vector_lo[0], &vector_hi[0]  );

#endif

  return(vcount);
}    

/* Method to check and reset the value for vector_chunk depending on the currently compiled
   in max memory parameter
*/

int PsociGAbasis::resetVectorChunk( int maxmempercore )
{
  int t_nchunk = 0;

if ( ! l_scatter_products ) {

/* a this point we could change the value of l_nxcheck to fit memory constraints */

static int orig_chunk = local_hamiltonian->fetchNchunk();

long estmem = local_hamiltonian->matrixVectorProductsMemoryPerCore();
long estmemMB = estmem/1000000.0;

if ( GA::nodeid() == 0 ) cout << " Total memory per core for M*v is (MB)" << estmem/1000000.0 << endl;
if ( GA::nodeid() == 0 ) cout << " Memory estimate only good for NEW CHUNK J method" << endl;

   const int maxmem=maxmempercore;

   if ( orig_chunk == 1 && estmemMB > maxmem ) {
      int frac = 1 + (estmemMB/maxmem);
      t_nchunk = orig_chunk * frac;
      local_hamiltonian->resetNchunk( t_nchunk );
      if ( GA::nodeid() == 0 ) cout << GA::nodeid() << " New chunk value is " << t_nchunk << "based on current maxmempercore of " << maxmem << endl;
      if ( GA::nodeid() == 0 ) cout << GA::nodeid() << " New memory reading is " << local_hamiltonian->matrixVectorProductsMemoryPerCore()/1000000.0 << endl;
   } else {
      if ( GA::nodeid() == 0 ) cout << "existing chunk is good or was overridden"<< orig_chunk << endl;
   }
  } // scatter

  return(local_hamiltonian->fetchNchunk());
}    


/** Take the current rnorms and create a single new vector adding it to the subspace.
    Orthonormalization is performed 
    cout << "CURRENT ALL VECTORS " << endl;
*/

// Add local_ga and set flag to compute new local_Hga
// Shift must be eval(curroot) in the current subspace

// Collective
int PsociGAbasis::preconditionAndNewVectorIncore( int curroot, double shift, GA::GlobalArray * g_rnorms, GA::GlobalArray * g_diag )
{
 
#ifdef DETAILEDCHECK
  char * Message = "preconditionAndNewVectorIncore:rnorms";
  g_rnorms->checkHandle( Message );
  
  char * Message2 = "preconditionAndNewVectorIncore:diag";
  g_diag->checkHandle( Message2 );
#endif
  
  // Fetch current root to construct new vector.
  int lo[2], hi[2];
  int mylo[2], myhi[2];
  
  if ( curroot > global_cols || curroot < 0 ) {
    cerr << "curroot specified is too big or too small " << endl;
    GA::error("preconditionAndNewVectorIncore: curroot specified is too big or too small ", curroot );
  }
  
#ifdef DETAILEDCHECK
  if (  curroot > CurrentNumSoughtRoots() ) {
    GA::error( "curroot > local_maxroots", curroot );
  }
#endif

  lo[0] = 0;
  hi[0] = global_rows - 1;
  lo[1] = curroot - 1;
  hi[1] = curroot - 1;
  
  //g_scratch
  mylo[0] = 0;
  myhi[0] = global_rows - 1;
  mylo[1] = 0;
  myhi[1] = 0;
  
  //fetch actual root, gets its rnorm and form a new vector for appending`
  
  //g_scratch->printDistribution(); // single update vector
  //local_ga->printDistribution(); // current set of vectors
  
  //g_scratch->zero(); //May not be required
  g_scratch->copyPatch('N', g_rnorms, lo, hi, mylo, myhi );
  
  precondition( shift, g_scratch, g_diag ); // Update scratch with (shifted) eigenvalues
  
  return(1); // change this to something useful
  
  // No need to orthonormalize at this point since the "append" method will already do this for us
}

// Collective: g_scratch is changed on exit
void PsociGAbasis::precondition( double shift, GA::GlobalArray * g_scratch, GA::GlobalArray * g_diag )
{
  /* Now do the arithmetic to construct a new vector
     Step One: the Davidson preconditioning. Could also consider
     doing the Subspace Projected Approximate Matrix (SPAM) procedure
     if we had a suitable approximation for H.
     
     We should move away from this anyway and either do SPAM or 
     Jacobi-Davidson
     
     Data-parallel approach
  */
  
  // Take advantage of the fact that ihi-ilo indexes are contiguous: Is this guarenteed by GA ?
  // It will not break but performance could be very bad if they are not contiguous
  
  int lo[2], hi[2];
  g_scratch->distribution( GA::nodeid(), lo, hi );

  if (lo[0] < 0 || lo[1] < 0 ) { 
     cerr << "precondition No scratch  data on me " << GA::nodeid() << endl;
      GA::Terminate();
   } 
  
  int width = hi[0] - lo[0] + 1;

  vector<double> vec( width, 0);
  vector<double> diag( width ,0);
  
  int n1=1;
  g_scratch->get( lo, hi, &vec[0], &n1);
  g_diag->get( lo, hi, &diag[0], &n1);
  
  for(int i=0; i< width; ++i ) 
    {
    if ( abs(diag[i]-shift) > DAVIDSON_MIN  ) {
    vec[i] /= (diag[i]-shift);
    }
      //vec[i] = (abs(diag[i]-shift) <= DAVIDSON_MIN ) ? vec[i]: vec[i] / (diag[i]-shift) ; 
      //cout << "vec shift and diag are " << vec[i] << " " << shift << " " << diag[i] << endl;
    }
  g_scratch->put( lo, hi, &vec[0], &n1);
  
#if 0
    cout << "DATA AFTER PRECONDITION" << endl;
    vector<double> rvec( global_rows, 0.0 );
    int ld=1;
    lo[0]=0;
    lo[1]=0;
    hi[0]=global_rows-1;
    hi[1]=0;
    g_scratch->get( lo, hi, &rvec[0], &ld);
    cout << "DUMP PRECONDITION BEFORE ORTHO NORM COEFS for 1" << endl;
    for(int k=0; k< global_rows; ++k ){
    cout << setprecision(15) << rvec[k] << endl;
    }
#endif

  GA::SERVICES.sync();
  //GA::sync(); // ensure everyone finished preconditioning their piece of vec
}
void PsociGAbasis::setMatrixVectorChunk( int nchunk )
{
     local_hamiltonian->setVectorChunk( nchunk );
     local_hamiltonian->printVectorChunk();
}

void PsociGAbasis::prefetchCIMAT()
{
     local_hamiltonian->fetchLocalHamiltonianData();
}

void PsociGAbasis::prefetchListCIMAT()
{
     local_hamiltonian->fetchListHamiltonianData();
}

void PsociGAbasis::generateListCIMAT()
{
     local_hamiltonian->generateListHamiltonianData();
}

void PsociGAbasis::destroyGACIMAT()
{
     local_hamiltonian->destroyGAspace();
}

void PsociGAbasis::destroyLocalCIMAT()
{
     local_hamiltonian->destroyLocalHamiltonianData();
}

void PsociGAbasis::generateStaggeredChunkList(int l_vectorchunk )
{
  local_hamiltonian->generateStaggeredChunkList( GA::nodeid(), l_vectorchunk );
}

void PsociGAbasis::destroyStaggeredChunkList()
{
  local_hamiltonian->destroyStaggeredChunkList();
}

// Choose one or the other
void PsociGAbasis::allocateScratchVectorSpace()
{
     local_hamiltonian->allocateLocalScratchMatrixVector();
}

void PsociGAbasis::allocateScratchVectorSpaceScatterAccumulate()
{
     local_hamiltonian->allocateLocalScratchMatrixVector();
     local_hamiltonian->resizeBufferScatterMethod();
}

// Can be used for either choice
void PsociGAbasis::destroyScratchVectorSpace()
{
     local_hamiltonian->destroyLocalScratchMatrixVector();
}

void PsociGAbasis::allocateScratchCIMATSpace()
{
     local_hamiltonian->allocateLocalCIMATarray();
}


 
  









      
    




