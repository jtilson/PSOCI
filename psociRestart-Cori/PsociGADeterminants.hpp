/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociGADeterminants:

 Classes: 

 Description: 

 History:


**************************************************************************************/
/**
 *   @file PsociDeterminants.hpp
 *
 */

#ifndef PSOCIGADETERMINANTS_H
#define PSOCIGADETERMINANTS_H

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "ga++.h"

#include "PsociDeterminants.hpp"

using namespace std;

/* We could use SHORTs or smaller for many of these entries BUT we eventualy want to increase substantially the number
   of basis functions and orbs. So let's not worry about space yet. Besides, this allocation is transient. I.e., 
   it will copied into GA and deleted - PRIOR to hamiltonian construction
*/

/* Really do not need full INTS for most of the data in this array. Nearly all data are 0,1,2 in value. Further
   The orbital values will not likely go higher than 32 K.
*/

// DANGER must make consistent with PsociDeterminants.hpp
#ifdef USESHORT
typedef short COMPRESS_SIZE;
#else  
typedef int COMPRESS_SIZE; //GA doesn't support shorts at this time
#endif

// BUFFER gets stored into a GA in the form D(1:maxlength, i,i) i == ith spatial function
// Current form of data for assembly and passage to a GA determinant array

/* The total actual size/order of this array is 
   all entries of type COMPRESS_SIZE

      mxopn
      nbf
      nelec
      ksym;
      jsym;
      nsefi;
      num_singly;
      num_doubly;
      ndeti;
      index;   // Absolute configuration index useful for accessing GA space
      doubly(1-numdoubly)
      singly(1-numsingly)
      spin_projection(1-ndeti)
      phase(1-ndeti)
      sumsa(1-ndeti)
      sumsb(1-ndeti)
repeat ndeti times:
      nalpha
      alpha_spins(1-nalpha)
repeat ndeti times:
      nbeta
      beta_spins(1-nbeta)

//Total sizes in order

   10 * sizeof( COMPRESS_SIZE ) +
   1 * nclosed * sizeof( COMPRESS_SIZE ) + //doubly
   1 * nopen * sizeof( COMPRESS_SIZE ) +   /singly
   4 * ndeti * sizeof( COMPRESS_SIZE ) + // spin,phase,sumsa,sumsb
   ndeti * (
   1 * sizeof( COMPRESS_SIZE ) +
   1 * nalpha * sizeof( COMPRESS_SIZE ) +
    )
   ndeti * (
   1 * sizeof( COMPRESS_SIZE ) +
   1 * nbeta * sizeof( COMPRESS_SIZE ) +
    )

   nbeta = nelec - nalpha

   mxopn is specified at the start of the calculation. Maximum num open shells
   mxocc is known = ( nelec + mxopn ) / 2. Maximum number of occupied orbitals
   mxclosed = nelec / 2; Maximum possible number cloised shell orbs
   mxdeti determined from the problem. Maximum number of determinants per spatial
   mxalpha <= mxocc 
   mxbeta <= mxocc
  
   The maximum size needed for creation of a properly sized GA then becomes:

   10 * sizeof( COMPRESS_SIZE ) +
   1 * mxclosed * sizeof( COMPRESS_SIZE ) +
   1 * mxopn  sizeof( COMPRESS_SIZE ) +
   4 * mxdeti * sizeof( COMPRESS_SIZE ) +
   2 * mxdeti * sizeof( COMPRESS_SIZE ) +          //nalpha/nbeta
   1 * mxdeti * mxocc * sizeof( COMPRESS_SIZE ) +
   1 * mxdeti * mxocc * sizeof( COMPRESS_SIZE )

   =  ( 2*mxdeti*mxocc + 6*mxdeti + mxclosed + mxopn + 10 ) * sizeof( COMPRESS_SIZE )

   mxopn    -> from the configs input file
   mxdeti   -> from PsociDeterminants.global_maxdet_perspatial
   mxclosed -> nelec / 2
   mxocc    -> ( nelec + mxopn ) / 2

*/

const int MAXWIDTH_CONFIGS=100000; // only useful to singlenode configs processing

class PsociGADeterminants {
  
private:
  GA::GlobalArray * g_det;
  GA::GlobalArray * g_det_num;
  GA::GlobalArray * g_nsef;
  GA::GlobalArray * g_orbmap;

  PsociConfigs * local_config;
  PsociDeterminants * local_det;

  vector<char> replicated_phases; // should be maxspatials long

  int l_lrow; // Local determinant patch coordinates
  int l_hrow;
  int l_lcol;
  int l_hcol;
  int l_maxlength; //Max length of a column in g_det
  int l_maxSpatials; // Need this if we tearDownDets after upload to GA
  int l_maxsef; 
  int l_maxdet;
  int l_nbfn;
  int l_nelec;
  int l_ksym;
  int l_maxdetperconf;
  int l_maxsefperconf;
  int l_mxopn;

// New tricks for speed in the download section of the code

  vector<int> l_doubly;
  vector<int> l_singly;
  vector<int> l_spin_projection;
  vector<int> l_phases;
  vector<int> l_sumsa;
  vector<int> l_sumsb;

  vector<vector<int> > alpha_spins;
  vector<vector<int> > beta_spins;
  vector<int> aspins;
  vector<int> bspins;

  vector<COMPRESS_SIZE> buffer;

  //int packDeterminants( JOUTFG & joutfg );
  int uploadBufferToGA( int index, int tstsize, vector<COMPRESS_SIZE> & buffer);
  int downloadBufferFromGA( int index, int mxlength, int &obslength, vector<COMPRESS_SIZE> &buffer );
  int downloadBufferFromGA2D( int index1, int index2, int maxsize, int & obs_size, vector<COMPRESS_SIZE> & buffer);
  int downloadBufferFromGAPreallocate( int index, int mxlength, int &obslength, vector<COMPRESS_SIZE> &buffer );
  int downloadBufferFromGA2DPreallocate( int index1, int index2, int maxsize, int & obs_size, vector<COMPRESS_SIZE> & buffer);
  int computeLength( int deti, int occ, int closed, int opn );
  void setParameters( PsociDeterminants * deters);


public:
  PsociGADeterminants( PsociConfigs * configs, GA::GlobalArray * global_det, GA::GlobalArray * global_det_num , GA::GlobalArray * global_nsef, GA::GlobalArray * global_orbmap  );
  PsociGADeterminants(  GA::GlobalArray * global_det, GA::GlobalArray * global_det_num, GA::GlobalArray * global_sefi_start, GA::GlobalArray * global_orbmap   );
  void createDistributedDeterminants( PsociDeterminants * deters );
  void initDistributedDeterminants( PsociDeterminants * deters );
  void createDistributedDeterminants(  PsociDeterminants * deters, pair<int,double> & info );
  GA::GlobalArray * fetchDeterminantHandle();
  GA::GlobalArray * fetchDeterminantNumberHandle();
  void printDistribution();
  void fetchNumDimensions( int & rows, int & cols );
  void printDimensions();
  void assembleDimensions();
  void assembleAndBrdcstNsef( vector<int> & l_nsef );
  int fetchMaxLength();
  int packAndLoadDeterminantData();
  int fetchAndUnpackDeterminantData( int index, JOUTFG & fetch );
  int fetchAndUnpackDeterminantData( int index, pair<int,double> & info, JOUTFG & fetch  );
  int fetchAndUnpackDeterminantData2D( int index1, int index2, vector<JOUTFG> & joutfg  );
  int fetchAndUnpackDeterminantDataPreallocate( int list, pair<int,double> & info, vector<COMPRESS_SIZE> & buffer, JOUTFG & joutfg  );
  int fetchAndUnpackDeterminantDataPreallocate( int list, vector<COMPRESS_SIZE> & buffer, JOUTFG & joutfg  );
  int fetchAndUnpackDeterminantData2DPreallocate( int index1, int index2, vector<COMPRESS_SIZE> & buffer, vector<JOUTFG> & toutfg  );


  int fetchNumElectrons();
  int fetchNumBasisFunctions();
  int fetchGlobalMaxSefPerSpatial();
  int fetchGlobalMaxDetPerSpatial();
  int fetchGlobalMaxSefs();
  int fetchGlobalSymmetry();
  int fetchMaxSpatials();
  void localConfRange( int & lo, int & hi );
  void uploadNsefToGA( int index, int nsefi );
  void printSpatialsData( JOUTFG & ioutfg );
  long generateDistributedDeterminants();
  long generateDistributedDeterminants(int crap);
  long generateDistributedDeterminantsSmallMemory( PsociConfigs & configs );
  int constructConfigurationData( JOUTFG & joutfg, string & configs );
  void preallocateTempJoutfgObject( JOUTFG & toutfg );
  void destroyJoutfgGA();
  PsociDeterminants * fetchInternalDeters();
  void destroyLocalDeterminantsData();
  long generateDistributedDeterminantsSmallMemoryDistributed( PsociConfigs & configs, vector<int> & list );
  void uploadSubsetDeterminants();
 // long generateDistributedDeterminantsSmallMemoryGA( PsociConfigs & configs );
 // long generateDistributedDeterminantsSmallMemoryGA(PsociConfigs & configs, pair<int,double> & info );
 long generateDistributedDeterminantsSmallMemoryGA();
 long generateDistributedDeterminantsSmallMemoryGA(pair<int,double> & info );
  int uploadOrbMaptoGA( int index,  vector<COMPRESS_SIZE> orbmap );
  void fetchOrbMap( int index, vector<COMPRESS_SIZE> & orbmap );
  void fetchOrbMap( int index, vector<vector<COMPRESS_SIZE> > & orbmap2d );
  int compareOrbMap( vector<COMPRESS_SIZE> & imap, vector<COMPRESS_SIZE> &jmap );
  int compareOrbMapMVP( vector<COMPRESS_SIZE> & imap, vector<COMPRESS_SIZE> &jmap );
  int compareOrbMapMVP( COMPRESS_SIZE * imap, COMPRESS_SIZE * jmap );

  int compareOrbMapCIDEN( COMPRESS_SIZE * imap, COMPRESS_SIZE * jmap );

  int compareAllOrbMap( int index, vector<COMPRESS_SIZE> & imap, vector<int> & list);
  int compareAllOrbMapMatrixVectorProducts( int index, vector<COMPRESS_SIZE> & imap, vector<int> & list);
  int compareAllOrbMapMatrixVectorProducts( int index, vector<COMPRESS_SIZE> & imap, vector<int> & jphase, vector<int> & list );

  int compareAllOrbMapCIDEN( int index, vector<COMPRESS_SIZE> & imap, vector<int> & list);
  int compareOrbMap( COMPRESS_SIZE * imap, COMPRESS_SIZE * jmap );
  int fetchAllCompressedOrbMap( vector<vector<pair<short,short> > > & ijpairlist );
  int orbmapRemoveZeros( vector<COMPRESS_SIZE> & orbmap, vector<pair<short,short> > & densepairlist );
  int compareAllOrbMapRemoveZeros( int index, vector<vector<pair<short,short> > > & ijpairlist, vector<int> & list );
  int compareOrbMapRemoveZerosCIDEN( int iindex, int jindex, vector<vector<pair<short,short> > > & ijpairs );
  int compareOrbMapRemoveZeros( int iindex, int jindex, vector<vector<pair<short,short> > > & ijpairs );
  int compareAllOrbMapRemoveZerosCIDEN( int index, vector<vector<pair<short,short> > > & ijpairlist, vector<int> & list );
  void computeConfigsGA( int index, vector<COMPRESS_SIZE> & orbmap, JOUTFG & toutfg );
  int compareLocalJSubsetOrbMap(int numJs, int index, vector<COMPRESS_SIZE> & imap, vector<int> & jindexlist, vector<COMPRESS_SIZE> & jsubsetlist, vector<int> & list );
  int compareLocalJSubsetOrbMap(int numJs, int index, vector<char> * rep_phase, vector<COMPRESS_SIZE> & imap, vector<int> & jindexlist, vector<COMPRESS_SIZE> & jsubsetlist, vector<int> & list );
  int compareLocalJSubsetOrbMapCIDEN(int numJs, int index, vector<COMPRESS_SIZE> & imap, vector<int> & jindexlist, vector<COMPRESS_SIZE> & jsubsetlist, vector<int> & list );
  int compareLocalJSubsetOrbMapCIDEN(int numJs, int index, vector<char> * rep_phase, vector<COMPRESS_SIZE> & imap, vector<int> &  jindexlist, vector<COMPRESS_SIZE> & jsubsetlist, vector<int> & list );

  void fetchOrbMap( int index, COMPRESS_SIZE * orbmap );
  vector<char> * fetchReplicatedPhasesHandle();


};

#endif //PSOCIGADETERMINANTS_H
