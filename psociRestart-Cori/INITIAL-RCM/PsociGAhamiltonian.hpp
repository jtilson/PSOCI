/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociGAhamiltonian:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociGAhamiltonian.hpp
 *
 */

#ifndef PSOCIGAHAMILTONIAN_H
#define PSOCIGAHAMILTONIAN_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <list>

#include <ga++.h>

#include "PsociDeterminants.hpp"
#include "PsociGADeterminants.hpp"
#include "PsociHamiltonian.hpp"
#include "PsociProcessMOcoefs.hpp"


extern "C" {
#include "mpi.h"
#ifdef USEMPI
#include "ga-mpi.h"
#endif
}

using namespace std;

const int MAX_CI_VECTOR_DIMS = 2; // Should equal MAX_GA_VECTOR_DIMS PsociGAbasis

const int NXTVAL_CHUNK_SIZE = 1;

const double  MIN_DENOMINATOR = 0.001 ; //Denominator check in precondition

const double MIN_DOUBLE_EVAL_PRINT=0.000000001;

const int MIN_HVEC_TOKEEP=0.0; // deflation method 

const int default_vector_chunk = 1;

const int MAX_AGGREGATE = 1; // Dont get too big:  maximum width of prefetched cimat/icicol/number also used for bufV2 and friends
const int MINIMAP_MAX = 2000; // 
const int MAX_DEFLATE_ELEMENTS=1;

//bool compare_evalues( pair<double,vector<double> > left, pair<double,vector<double> > right );


struct VECTRIPLE {
   int row;
   int col;
   double value;
};

struct  vector_vectortriplets
{
bool operator() ( const VECTRIPLE & left, const VECTRIPLE & right ) { return left.row < right.row; }
};

bool compare_vectortriplets( const VECTRIPLE & left, const VECTRIPLE & right );
bool compare_vectorpairs( const pair<int,double> & left, const pair<int,double> & right );


class PsociGAhamiltonian {
  

private:

  int l_nbf;
  int l_nelec;
  int l_numconfs;
  int l_maxspatials;
  int l_maxsparse_supp;
  int l_maxwidth;

  int l_nxtval_chunk;

  int l_nchunk;

  int l_granularity; // Parameter to indicate degree of bunching up of GA calls for determinant data
  int l_numJchunksHbuild;// 

// Hold local patch of H on the core
  vector<vector<int> > local_icicol;
  vector<vector<double> > local_cimat;
  vector<int> local_number;
  vector<int> cimat_list;

// Hold local scratch for bufV, scatter/gather arrays

// These are (potentially) allocated in allocateLocalScratchMatrixVector()
    vector<double> bufV; // #ifndef NEWACCUMULATEDSCATTER
    vector<double> dataarray;

/* Only used foir code directly calling nga_scatter_acc. Calling lowlevel pnga_scatter_acc
   doesn;t require this
*/ 
    int **indexarray2d; // experimental

    vector<Integer> indexarrayInteger; // choose one or the other at M*v time
    vector<int> indexarray;

    vector<double> vectors;


// These are potentially used in matrixVectorProductsSplitChunkScatterAccumulate
// They are resived in a public method
// ifdef NEWACCUMULATEDSCATTER

    vector<vector<double> > bufV2;
    vector<vector<int> > bufVIndex;
    vector<int> bufVnumber;

// All the above are cleared in destroyLocalScratchMatrixVector()
/*
   dataarray.clear();
   indexarray.clear();
   vectors.clear();
   bufV.clear();
   bufV2.clear();
   bufVIndex.clear();
   bufVnumber();
*/

 
// Hold potential staggers list
   vector<int> my_chunklist;

/* The following dims are needed in matrix-vector products routine. No need to 
   constantly reread it.
*/
  int local_dims[ 2 ]; 
  int local_cilo[2], local_cihi[2];


  GA::GlobalArray * l_g_cimat; // declare global det data
  GA::GlobalArray * l_g_icicol;
  GA::GlobalArray * l_g_number;
  GA::GlobalArray * l_g_diag_sef;

/* New bit for building a vector mark */

  GA::GlobalArray * l_g_mask;
  GA::GlobalArray * l_g_mask_number;


/* An experimental approach to decreasing aggregate GA memory being used */
// These are not necessarily used

  GA::GlobalArray * l_g_cimat_supp; 
  GA::GlobalArray * l_g_icicol_supp;
  GA::GlobalArray * l_g_number_supp;

  GA::GlobalArray * g_suppval;


  PsociHamiltonian  * l_gaHamiltonian;
  PsociGADeterminants * l_deters;

  PsociGADeterminants * l_deters2; // new for transition moments

/* New bits simply for getting H matrix data
*/
  fstream ohamiltonianfile;
// End extra bits

  vector<JOUTFG> l_confs; // Only used for the debugging method: localHamiltonian 

  vector<int> l_nsef;

/* THe next stuff is for the subspace testers and will soon disappear */

  long test;
  int global_rows, global_cols; // Reports actual size of matrix space
  GA::GlobalArray * local_ga;

  GA::GlobalArray * g_nxtval;

  vector<double> sefsef;


  int l_maxsparse;
  int l_maxsef;
  
  void createGAhamiltonian();
  void createGAhamiltonian2Array();
  void createGAhamiltonianDiagOnly();
  int uploadCIMATtoGA( int isef , vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number );
  int uploadCIMATtoGAAcc( int isef , vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number );
  int uploadCIMATtoGASupplemental( int isef , vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number );
  int uploadCIMATtoGASupplementalAcc( int isef , vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number );

  long numberNonzeroElements( long l_number );
  long nxtask( int nproc ); 
  long nxtaskAlt( int nproc );
  void initNxtask();
  void destroyNxtask();
  void aggregateStatsData(int l_maxsef, vector<vector<int> > & matrix );
  void zeroVector( vector<double> & vector );
  void zeroVector( vector<double> & vector, int number );
  void zeroVector( vector<int> & vector );
  void zeroVector( vector<int> & vector, int number );
  void zeroMatrix( vector<vector<double> >& matrix );
  void zeroMatrix( vector<vector<int> >& matrix );
  void zeroMatrix( vector<vector<double> >& matrix, int number );
  void zeroMatrix( vector<vector<int> >& matrix, int number );

   
  int fetchFullCIvector(int index, GA::GlobalArray * g_v, vector<double> & vect );
  int fetchFullCIvectorBRD(int index, GA::GlobalArray * g_v, vector<double> & vect );
  int fetchPartialCIvector(const int iroot, const int jstart, const int nsefj, GA::GlobalArray * g_v, vector<double> & vec_j );

  inline int twoDindex( int m, int n );
  void outputPacked2Dmatrix(bool square, int nbf, vector<double> & array, string title );
  void outputPacked1Darray(int nelem, vector<double> & array, string title );
  void outputPrintDensityAObasis( int nsym, vector<int> & nsize, vector<int> & nstart, vector<double> & density, string title );
  void replicateCIdensity(bool square, int nbf, int nroots, vector<vector<double> > & array );

  void transformNOstoAObasis( int isym, vector<double> & no_occupations, vector<double> & dens_vector_sym, vector<double> & dens_vectorAO_sym, PsociProcessMOcoefs & coefs );
  void sortNOinMObasis(int isym, vector<double> & evals, vector<double> & vectors );
  int fetchListWithinRange( int start, int end, vector<int> & icicol, int icinum, vector<int> & list );
/* not reallu local change namne */
  void fetchLocalHblock( int i, vector<double> & cimat, vector<int> & icicol, int & number );
  void fetchAllLocalHblock( int clow, int chi, vector<vector<double> > & local_cimat, vector<vector<int> > & local_icicol, vector<int> & local_number );
  void fetchListHblock( vector<int> & list, vector<vector<double> > & local_cimat, vector<vector<int> > & local_icicol, vector<int> & local_number );
  int gatherVectorSet(GA::GlobalArray * g_v, const int & ivec, vector<int> & list, vector<double> & vectors,  vector<double> & data, vector<int> & indexarray );

  void fetchAndReplicateVectorChunk( int ichunk, int ivec, int width, int * v_lo, int * v_hi, vector<double> & chunk, GA::GlobalArray * g_v );

  void accumulateAndPutProductChunk( int jchunk, int ivec, int width, int * v_lo, int * v_hi, vector<double> & bufV, GA::GlobalArray * g_Hv );

  void fetchAndReplicateVectorChunkMPI(int ichunk, int ivec, int width, int * v_lo, int * v_hi, vector<double> & chunk, GA::GlobalArray * g_v );

  void replicateAndPutHV( int ivec, int maxsefs, int * v_lo, int * v_hi, vector<double> & bufv, GA::GlobalArray * g_Hv );

public:

//  PsociGAhamiltonian( void );
  PsociGAhamiltonian( string & dirstring, int maxsparse, int maxsef,int maxsparseBig, int maxwidth, PsociGADeterminants * deters, PsociIntegrals * ints);
  PsociGAhamiltonian(int maxsparse, int maxsef,int maxsparseBig, int maxwidth,PsociGADeterminants * deters,PsociIntegrals * ints);
  PsociGAhamiltonian(int maxsparse, int maxsef,PsociGADeterminants * deters,PsociIntegrals * ints);
  PsociGAhamiltonian(int maxsparse, int maxsef, int supp_maxsparse, int maxwidth,  PsociGADeterminants *deters_det, PsociGADeterminants * deters_det2, PsociIntegrals * mos );

  ~PsociGAhamiltonian();
  void printKsym();
  void printPerSpatialValues();
  int setLocalConfigs( vector<JOUTFG> & confs );
/*
  int constructLocalHamiltonian();
  long constructGlobalHamiltonian();
  long constructGlobalHamiltonianExp();
  long constructGlobalHamiltonianExp( pair<int,double> & info );
  long constructGlobalHamiltonianDirect( pair<int,double> & info );
  long constructGlobalHamiltonianDirect();
  long constructGlobalHamiltonianExpDirect( pair<int,double> & info );
  long constructGlobalHamiltonianExpDirect();
  long constructGlobalHamiltonianExpGood( );
  long constructGlobalHamiltonianExpAll();
  long constructGlobalHamiltonianExpAll( pair<int,double> & info );
  long constructGlobalHamiltonianExpAlternative();
  long constructGlobalHamiltonianExpAlternative( pair<int,double> & info );
  long constructGlobalHamiltonian2WAY( pair<int,double> & info );
  long constructGlobalHamiltonian2WAY();
  long constructGlobalHamiltonian2WAYDirect( pair<int,double> & info );
  long constructGlobalHamiltonian2WAYDirect();
  long constructGlobalHamiltonian( pair<int,double> & info );
  long constructGlobalHamiltonianSparseStats();
  long constructGlobalTransitionHamiltonian2WAY(pair<int,double> & info);
  long constructGlobalTransitionHamiltonian2WAY();
*/

  //long constructGlobalHamiltonianForCIDEN2WAY( int numRoots,  GA::GlobalArray * g_v, vector<vector<double> > & density);
  //long constructGlobalHamiltonianForCIDEN2WAY(  int numRoots,  GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & info );
  
 // long constructGlobalHamiltonianForCIDEN2WAY();
  //long constructGlobalHamiltonianForCIDEN2WAY(pair<int,double> & info );
  
// New stuff


  long constructGlobalTransitionHamiltonian2WAY(pair<int,double> & info);
  long constructGlobalTransitionHamiltonian2WAY();
  long constructGlobalHamiltonian2WAY( pair<int,double> & info );
  long constructGlobalHamiltonian2WAY();
   long constructGlobalHamiltonianAltOrbMapDirect(pair<int,double> & info );
   long constructGlobalHamiltonianAltOrbMapDirect();
   long constructCIdensityAltOrbMapDirect( pair<int,double> & info );
   long constructCIdensityAltOrbMapDirect();

// end new stuff

  void flushH();

  int constructMASKGAAcc(int * v_lo, int * v_hi, GA::GlobalArray * g_v);
  int fetchMaskList( int my_rank, vector<int> & mask_list );
  int matrixVectorProductsIncore(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );
  int matrixVectorProductsIncore(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims , int * v_lo, int * v_hi ,pair<int,double> & info );
  int matrixVectorProductsIncoreGOP(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );
  int matrixVectorProductsIncoreGOP(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims , int * v_lo, int * v_hi ,pair<int,double> & info );

  int matrixVectorProductsIncoreGOPSplit(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

 int matrixVectorProductsIncoreGOPSplit(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info );

 int matrixVectorProductsIncoreGOPSplitChunk(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info );

 int matrixVectorProductsIncoreGOPSplitChunk(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

 int matrixVectorProductsIncoreGOPSplitChunkScatter(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info );

 int matrixVectorProductsIncoreGOPSplitChunkScatter(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

 int matrixVectorProductsIncoreGOPSplitGatherScatter(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );
 int matrixVectorProductsIncoreGOPSplitGatherScatter(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info );

 int matrixVectorProductsIncoreGatherScatter(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );
 int matrixVectorProductsIncoreGatherScatter(int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims , int * v_lo, int * v_hi ,pair<int,double> & info );


 long generateIEXJEXforCIDEN(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iex, vector<int> & jex );


/* More subspace stuff */

  int CurrentNumRows();
  int CurrentNumCols();
  GA::GlobalArray * fetchHandle();
  int CurrentNumSoughtRoots();
  int testBuild();
  int computeNsefEntryPoints( vector<int> & nsef ); 
  int packSefsef( JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef, vector<vector<double> > & cimat, vector<vector<int> > & icicol, vector<int> & number );
  int packSefsef( JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef, vector<vector<double> > & cimat, vector<vector<int> > & icicol, vector<int> & number,  
		  vector<pair<long,long> > & stats, vector<vector<int> > & matrix );
  void printNsefEntryPoints();
  void hamiltonianDistribution();
  void dumpCimat();
  int fetchListOfDiags( int lo, int hi, vector<double> & diags );
  void printDiags();
  int fetchNumberSefs();
  int fetchDiags( vector<double> & diags );
  int precondition( double shift, GA::GlobalArray * g_v, vector<double> & diags );
  
  
  void initSuppVal();
  void destroySuppVal();
  long nextSupplimentIrow();
  void SetGranularity( int width );
  void PrintGranularity();
  int fetchNxtvalChunk(void );
  void setNxtvalChunk( int nchunk );
  int fetchNBF();
  GA::GlobalArray * fetchCIMAT(); 
  GA::GlobalArray * fetchICICOL(); 
  GA::GlobalArray * fetchNUMBER(); 
  GA::GlobalArray * fetchDIAG_SEF();
  GA::GlobalArray * fetchCIMAT_supp();
  GA::GlobalArray * fetchICICOL_supp();
  GA::GlobalArray * fetchNUMBER_supp();
  void freeHamiltonian();
  int transformSEFSEFtoMO( int numRoots, JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef, GA::GlobalArray * g_v, vector<vector<double> > & density );
  int transformSEFSEFtoMOnoGET(JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef, vector<double> & vec_i, vector<double> & density );
  int transformSEFSEFtoMOGET(int iroot, JOUTFG & iconf, JOUTFG & jconf , vector<double> & sefsef, GA::GlobalArray * g_v , vector<double> & density );
  //long constructCIdensity( int numRoots, GA::GlobalArray *  g_v, vector<vector<double> > & density, pair<int,double> & CIDinfo  );
  //long constructCIdensityAlt( int numRoots, GA::GlobalArray *  g_v, vector<vector<double> > & density, pair<int,double> & CIDinfo  );
  //long constructCIdensity( int numRoots, GA::GlobalArray *  g_v, vector<vector<double> > & density);
  //long constructCIdensityAlt( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density );
  long constructCIdensityAltOrbMap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density );
  long constructCIdensityAltOrbMap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & CIDinfo  );
  int fetchSEFSEF( JOUTFG & fetchi, JOUTFG & fetchj, vector<double>  & sefsef );
  int fetchNELEC();
  long constructCIdensityAltOrbMapDirect( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density );
  long constructCIdensityAltOrbMapDirect( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & CIDinfo  );
  //long constructCIdensityAltOrbMap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density );
  //long constructCIdensityAltOrbMap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & CIDinfo  );

  long constructCIdensityAltOrbMapDirectChunkJmap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density );
  long constructCIdensityAltOrbMapDirectChunkJmap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & CIDinfo  );

  long constructCIdensityAltOrbMapChunkJmap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density );
  long constructCIdensityAltOrbMapChunkJmap( int numRoots, GA::GlobalArray * g_v, vector<vector<double> > & density, pair<int,double> & CIDinfo  );

  void setVectorChunk( int nchunk );
  void printVectorChunk();
  
   void fetchLocalHamiltonianData();
   void fetchListHamiltonianData();
   void fetchsubListHamiltonianData(vector<int> & list );
   void generateListHamiltonianData();

   void destroyLocalHamiltonianData();
  void generateStaggeredChunkList(int nodeid, int nchunk );
  void destroyStaggeredChunkList();
  int deflateMultimapToVector( multimap< int, double > & map, vector<double> & data, vector<int> & indexarray, const int & ivec );
//  int deflateMultimapToVector( multimap< pair<int,int>, double > & map, int deflate_index, vector<double> & data, vector<int> & indexarray );


  int deflateListToVector( list<VECTRIPLE> & vector_buf, vector<double> & data, vector<int> & indexarray);
  int deflateVectorToVector( vector<VECTRIPLE> & vector_buf, vector<double> & data, vector<int> & indexarray);
  void destroyGAspace();

int matrixVectorProductsIncoreGatherScatterNoVectorGet(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info ); 
int matrixVectorProductsIncoreGatherScatterNoVectorGet( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv,
                                                                 int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsIncoreGatherScatterNoVectorGetDecember( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv,
                                                                 int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsIncoreGatherScatterLocalVectorGetMarch( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv,
                                                                 int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsSplitChunkRepAccumulate(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info );
int matrixVectorProductsSplitChunkRepAccumulate(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsSplitLocalVectorRepAccumulate(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );


int matrixVectorProductsSplitChunkScatterAccumulate( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsSplitChunkScatterAccumulate(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info );

int matrixVectorProductsSplitChunkScatterAccumulateLowLevel(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi , pair<int,double> & info );

int matrixVectorProductsSplitChunkScatterAccumulateLowLevel(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulate(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirect(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsSplitChunkScatterAccumulateLowLevelAccumulateDirectFetch(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );
int matrixVectorProductsSplitAccumChunkGatherVector(int start, int end, GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

int matrixVectorProductsSplitChunkMultipleJchunks( int start, int end,  GA::GlobalArray * g_v, int * v_dims, GA::GlobalArray * g_Hv, int * hv_dims, int * v_lo, int * v_hi );

 void allocateLocalScratchMatrixVector();
 void destroyLocalScratchMatrixVector();
 void allocateLocalCIMATarray();
 long matrixVectorProductsMemoryPerCore();

 void resizeBufferScatterMethod();
 int packDataNoZeros(int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, vector<int> & indexarray, vector<double> & dataarray );
 int **packDataNoZeros(int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, int & actualsize, vector<double> & dataarray );
 int packDataNoZerosTransposed(int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, vector<int> & indexarray, vector<double> & dataarray );

 int packDataNoZeros1D(bool transpose, int ivec, int number,  vector<int>& index, vector<double>& bufV, vector<int> & indexarray, vector<double> & dataarray );

 void emptyIndexarray( int num, int**data );
 void LoopNputMethod(GA::GlobalArray * g_a, vector<int> & indexarray, vector<double> & values, GANbhdl * g_handle );


 void indexToInteger( vector<int> & indexarray, vector<Integer> & indexarrayInteger );

 int fetchJsubsetlist( int jchunk, vector<int> & jlistindex, vector<COMPRESS_SIZE> & jsubsetlist, int chunk_width, int l_maxspatials );
 long constructGlobalHamiltonianAltOrbMapDirectChunkJmap();
 long constructGlobalHamiltonianAltOrbMapDirectChunkJmap( pair<int,double> & info );
 long constructGlobalHamiltonianAltOrbMapChunkJmap();
 long constructGlobalHamiltonianAltOrbMapChunkJmap( pair<int,double> & info );

 long constructGlobalHamiltonianAltOrbMapDirectChunkJmapInner();
 long constructGlobalHamiltonianAltOrbMapDirectChunkJmapInner( pair<int,double> & info );

 void SetJChunk( int numchunks );
 void PrintJChunk();

 int packDataNoZerosAccumulate(bool transpose, int indexnewearray, int ivec, vector<int> & number,  vector<vector<int> >& index, vector<vector<double> >& buf, vector<int> & indexarray, vector<double> & dataarray );

 int packDataNoZerosAccumulateSingle(bool transpose, int indexnewarray, int ivec, int  number,  vector<int> & index, vector<double> & buf, vector<int> & indexarray, vector<double> & dataarray );

 void fetchPatchVector( int ivec, int ch_start, int ch_end, vector<double> & testvector, GA::GlobalArray * g_v );
 void fetchAllVector( int ivec, vector<double> & testvector, GA::GlobalArray * g_v );

 int matrixVectorUpdate(int iconf, int & indexnewarray, int nsefi,  vector<int> & bufVnumber, vector<vector<int> > & bufVIndex, vector<vector<double> > & bufV2,vector<vector<int> > & icicol, vector<vector<double> > & cimat, vector<int> & number, vector<double> & testVector );

 int matrixVectorUpdate(GA::GlobalArray * g_v ,int iroot, int iconf, double Vextra, int * v_lo, int * v_hi, int & indexnewarray, int nsefi,  vector<int> & bufVnumber, vector<vector<int> > & bufVIndex, vector<vector<double> > & bufV2,vector<vector<int> > & icicol, vector<vector<double> > & cimat, vector<int> & number, vector<double> & testVector );

 int matrixProductScatter( int iroot, int & indexnewarray, int & estimate_elements, vector<int> & bufVnumber, vector<vector<int> > & bufVIndex, vector<vector<double> > & bufV2, vector<int> & indexarray, vector<double> & dataarray, GA::GlobalArray * g_Hv );

 int uploadCIMATdiagElemsOnlytoGA( int isef, vector<vector<double> > &cimat, vector<vector<int> > &icicol, vector<int> & number );
 long computeOnlyDiagElements();
 int assemblePhasePossibles(int iconf, int iphase, vector<int> & jphases, vector<char> * replicated_phases );
 int assemblePhasePossiblesCIDEN( int iconf, int iphase, vector<int> & jphases, vector<char> * replicated_phases );

 void resetNchunk( int nchunk );
 int fetchNchunk(); 

};

#endif //PSOCIGAHAMILTONIAN_H
