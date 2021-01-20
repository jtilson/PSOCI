
/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: This method computes the matrix element sub-blocks given two JOUTFG objects

              Generally these are LOCAL calls operating on local data

 History:
          First release: Sept 1, 2011
          Add preproceessing flag DETAILEDHAMCHECK to turn on detailed checking
          mostly used for debugging and method modification

          Corrected the return values on several methods so that 
          proper screening on the blocks can happen

**************************************************************************************/

/**
 *   @file PsociHamiltonian.C
 *
 */

#include <cmath> 

#include "PsociTimer.hpp"
#include "PsociDeterminants.hpp"
#include "PsociHamiltonian.hpp"

using namespace std;

PsociHamiltonian::PsociHamiltonian(int nbf, int nelec )
{
  l_nbf = nbf;
  l_nelec = nelec;
  l_maxdet = 0;
  l_maxsef = 0; 
  l_ksym = 0;
}

void PsociHamiltonian::printKsym()
{
  if ( GA::nodeid() == 0 ) {
  cout << "Current global symmetry is " << l_ksym << endl;
  }
}

void PsociHamiltonian::specifyIntegrals( PsociIntegrals * ints )
{
  l_integrals = ints;
}

void PsociHamiltonian::setPerSpatialValues( int maxdet, int maxsef )
{
  l_maxdet = maxdet;
  l_maxsef = maxsef;
} 

void PsociHamiltonian::setKsym( int ksym )
{
  l_ksym = ksym; //Global symmetry
}

void PsociHamiltonian::printPerSpatialValues()
{
  if ( GA::nodeid() == 0 ) {
  cout << "Current l_perspatial_maxdet is " << l_maxdet << endl;
  cout << "Current l_perspatial_maxsef is " << l_maxsef << endl;
  }
}


//Timer/Accumulator wrapped version
int PsociHamiltonian::generateOffdiagBlocks(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iocc, vector<int> & jocc, vector<int> & iex, vector<int> & jex, vector<double> & vnt, 
					    vector<double> & detdet , pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int processed = generateOffdiagBlocks( iconf, jconf, iocc, jocc, iex, jex, vnt, detdet );
  info.second += psociTime() - timein; //Accumulate times
  info.first = g_rank;
  return( processed );
}

/* This method performs multiple steps. 1) The identification of the excited orbitals (formally called iex and jex)
   for the input configurations i and j. and the computations of the excitation number. The calling program determines
   which values should trigger a matrix block calculation but generally >2 means no interactions.  2) for all determinants
   the block of nij arrays are computed
*/
int PsociHamiltonian::generateOffdiagBlocks(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iocc, vector<int> & jocc, vector<int> & iex, vector<int> & jex, vector<double> & vnt, 
					    vector<double> & detdet  )
{
#ifdef DETAILEDHAMCHECK
  if ( iex.size() == 0 || jex.size() == 0 ) {
    cerr << "OffdiagBlocks: iex or jex are empty: non-physical scenario: iex is " << iex.size() << " jex is " << jex.size() << endl;
    GA::Terminate();
  }
#endif
  
  //If greater than 4 no calculation required. But this was previously checked

  int deltaOpn =  abs( iconf.num_singly - jconf.num_singly );
  
 // vector<pair<int,int> > nmark; // nmarka,nmarkb for each pair
  
 // Determine excitation number and excited orbital list
 // If number > 2 then entire block interaction is zero by symmetry
  
  //Start processing excitations: We could prune down the iex/jex construction above by uncommenting out the breaks.
  
  int kk = iex.size();
  
 // vector<double> vnt; // iconf-jconf summed integral contribution
  
  int do_sefcig = 0;
  if ( kk == 2 ) { 
//    vnt.resize(2);
    do_sefcig = processDoubleExcitation( iex, jex, vnt );
    // Process single excitations
  } else if ( kk == 1 ) { 
//    vnt.resize(1);
    do_sefcig = processSingleExcitation( iconf, jconf, iocc, jocc, iex, jex, vnt, detdet );
  } //end excitations  
  
  return( do_sefcig );
}


/* Compute i<>j double excitation contributions
*/
int PsociHamiltonian::processDoubleExcitation( vector<int> & iex, vector<int> & jex, vector<double> & vnt )
{
  vnt.clear();
  int have_processed = 1; // Innocent
  
  /*
    double excitation
    
    the excitations can be classified according to (1) the difference
    in the number of open shells (deltaOpn) and (2) the sum of the
    occupation numbers of the excited orbitals (jmarr) as follows:
    
    excitation type
    
     now = 0       jmarr = 4           ab = cd
                         = 8           a(2) = b(2)
     now = 2       jmarr = 5           a(2) = bc
                         = 6           a(2)b = acd
     now = 4       jmarr = 6           a(2)b(2) = abcd
     
  determine the integrals needed and the alignment
  
  */

// Might want to disable these checks for speed

#ifdef DETAILEDHAMCHECK
  if ( iex.size() == 0 || jex.size() == 0 ) {
    cerr << "DoubleExcitation: iex or jex are empty: iex is " << iex.size() << " jex is " << jex.size() << endl;
    GA::Terminate();
  }
  if ( iex.size() != 2 || jex.size() != 2 ) {
    cerr << "Non conforming iex or jex: Abort: sizes: iex = " << iex.size() << " jex is " << jex.size() << endl;
    GA::Terminate();
  }
  if ( l_integrals == NULL ) {
    cerr << "Integrals not specified: Aborting " << endl;
    GA::Terminate();
  }
#endif
  
  vnt.push_back( l_integrals->fetch2eIntegralValue( iex[0], jex[0], iex[1], jex[1] ) );
  vnt.push_back( l_integrals->fetch2eIntegralValue( iex[0], jex[1], iex[1], jex[0] ) );
  
  const double dZero = 0.0;
  if ( vnt[0] == dZero && vnt[1] == dZero ) have_processed = 0;
  return( have_processed );
}

//On return if Zero then no i,j contributions were required or integrals were screened to zero 

int PsociHamiltonian::processSingleExcitation( JOUTFG & iconf, JOUTFG & jconf, vector<int> & iocc, vector<int> & jocc, vector<int> & iex, vector<int> & jex, vector<double> & vnt , 
					       vector<double> & detdet )
{
  vnt.clear();
  const double dZero = 0.0;
  int do_sefcig = 1; // sefcig: yes unless told otherwise

/*
    single excitation.
    determine vnt(1) and the alignment
    
     vnt(1) = (iex1,h,jex1) + sum(core)(iex1,2jc-kc,jex1)
                            + sum(open)(iex1,ji,jex1)
                            - (iex1,jiex1,jex1)
*/

#ifdef DETAILEDHAMCHECK
  if ( iex.size() == 0 || jex.size() == 0 ) {
    cerr << "processSingleExcitation: iex or jex are empty: iex is " << iex.size() << " jex is " << jex.size() << endl;
    GA::Terminate();
  }
#endif
  
  int iext = iex[0];
  int jext = jex[0];
  
  double vnt1;

  int jsymi = iconf.jsym;
  int jsymj = jconf.jsym;
  
  // spin-orbit matrix element

  if ( jsymi != jsymj ) {
    vnt1 = l_integrals->fetch1eIntegralValue( iext, jext );

    double dZero = 0.0;
    if ( vnt1 != dZero ) {
    if ( jext > iext ) vnt1 = -vnt1;
    vnt1 = vnt1 / 2;

    vnt.push_back( vnt1 );

//TODO If here then we do not want to call sefcig  

      socig(iconf, jconf, iex, jex, vnt, detdet);
    }
    do_sefcig = 2; // Special assertion for this method

  } else {

  vnt1 = l_integrals->fetch1eIntegralValue( iext, jext ) - l_integrals->fetch2eIntegralValue( iext, jext, iext, iext );
    
   double buf=0;
   for(int i=0; i< l_nbf; ++i) {
     if ( iocc[i] != 0 ) {
       buf = l_integrals->fetch2eIntegralValue( iext, jext, i, i );
       vnt1 += buf;
       if ( iocc[i] == 2 ) {
	 vnt1 += buf; // Two equivelent contributions
       }
       if (iocc[i] == 2 && jocc[i] == 2 ) {
	 vnt1 -= l_integrals->fetch2eIntegralValue( iext, i, i, jext );
       }
     }
   }
    vnt.push_back( vnt1 );
    if (vnt1 == dZero ) do_sefcig=0;
  }
  //have_processed = 1;
  return( do_sefcig );
}

/* This is a convenience function and will probably be little used
   Ensures but arguments are the same, then processes the first 
*/
int PsociHamiltonian::generateDiagonalBlocks(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iocc, vector<double> & vnt )
{
#ifdef DETAILEDHAMCHECK
  if ( iconf.completed_state !=1 || jconf.completed_state !=1 ) {
    GA::error("generateDiagonalBlocks: One of both input confs are not in completed_state: Aborting ",1);
  } 
  if (iconf.index != jconf.index ) {
    GA::error("generateDiagonalBlocks: Only processes i == j terms: Aborting ",1);
  }
  if ( iocc.size() == 0 ) {
    GA::error("generateDiagonalBlocks: empty iocc ",1);
  }
#endif
  
  int processed = generateDiagonalBlocks( iconf, iocc, vnt );
  return( processed );
}

//Timer wrapped version
int PsociHamiltonian::generateDiagonalBlocks(JOUTFG & iconf, vector<int> & iocc, vector<double> & vnt, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int processed = generateDiagonalBlocks(iconf, iocc, vnt);
  info.second += psociTime() - timein; //Accumulate times
  info.first = g_rank;
  return( processed );
}

/* Actual matrix processing of selected diag terms 
 */
int PsociHamiltonian::generateDiagonalBlocks(JOUTFG & iconf, vector<int> & iocc, vector<double> & vnt )
{  
  int processed = processDiagExcitation( iconf, iocc, vnt ); 
  return( processed );
}

/* Sum integrals for the H(iconf,iconf) term
 */

int PsociHamiltonian::processDiagExcitation( JOUTFG & iconf, vector<int> & iocc, vector<double> & vnt )
{
  int have_processed = 0;
  int nopeni = iconf.num_singly;
  
  if ( nopeni > 1 ) {
    vnt.resize( 1 + nopeni * (nopeni-1) / 2 ); //( open shell combinations plus the zeroth term ) 
    
    /* 
       kints obtained from transpositions on the open shells of
       configuration i.
    */ 
    
    int iopnii, iopnjj;
    int mm=1; // Yep we skip the first entry
    for(int i=1; i< nopeni; ++i ) { 
      for(int j=0; j< i; ++j ) {
	
	iopnii = iconf.singly[i]; 
	iopnjj = iconf.singly[j];
	if ( iopnii > l_nbf || iopnjj > l_nbf ) {
	  cerr << "Fatal error: iopnii or iopnjj > nbf. nbf = " << l_nbf << "iopnii , iopnjj " << iopnii << " " << iopnjj << endl;
	  GA::Terminate();
	}
	vnt[ mm ] = ( l_integrals->fetch2eIntegralValue( iopnii, iopnjj, iopnjj, iopnii ) );
        ++mm;
      }
    }
  } else {
    vnt.resize(1); //Need to allocate space for other circumstances
  }

/*
     determine vnt(1)

     vnt(1) = sum(core+opens)*(ii,h,ii)*occnum(ii)
                  + sum(core)*sum(open)*(2jci-kci)
                  + sum(c1.le.c2)*(2jc1c2-kc1c2)
                  + sum(c1.lt.c2)*(2jc1c2-kc1c2)
                  + sum(i1.lt.i2)ji1i2
*/

//TODO grab integrals in batches for all i,j. But only for relevant i and j.

  double dTwo = 2.0;
  double vnt1 = 0.0;
  int ninj;
  double buf;
  
#ifdef DETAILEDHAMCHECK
  if ( iocc.size() != l_nbf ) {
    cerr << "Error ir iocc size. is " << iocc.size() << endl;
    GA::Terminate();
  }
#endif
  
  for(int ibf = 0; ibf < l_nbf; ++ibf ) {
    if( iocc[ibf] != 0 ) {
      if( iocc[ibf] == 1) {
	vnt1 += l_integrals->fetch1eIntegralValue( ibf, ibf );
      } else {
	vnt1 += dTwo * l_integrals->fetch1eIntegralValue( ibf, ibf ) + l_integrals->fetch2eIntegralValue( ibf,ibf,ibf,ibf );
      }
      if ( ibf != 0 ) { // May not be needed
	for(int jbf = 0; jbf < ibf; ++jbf ) {
          if ( iocc[jbf] != 0 ) {
	    ninj = iocc[ibf] * iocc[jbf];
	    if ( ninj == 1 ) {
	      vnt1 += l_integrals->fetch2eIntegralValue( ibf,ibf,jbf,jbf );
	    } else {
	      buf = dTwo * l_integrals->fetch2eIntegralValue( ibf,ibf,jbf,jbf );
	      buf -= l_integrals->fetch2eIntegralValue( ibf,jbf,ibf,jbf );
	      ( ninj == 2 ) ? vnt1 += buf: vnt1 += (dTwo * buf);
	    }
	  }
        }
      }
    }
  }
  vnt[0] = vnt1; 
  have_processed = 1;
  return( have_processed );
}


//If return <> 0 then process through sefcig.
//If return < 0 we have an untrapped error

/* On entry we presume this is called by the GA-SOCI routines which populate the CI GA
   This gets triggered if no valid DRA restart file is found ( or if user overrides )
   This particular routine computes all I/J double group terms  into the submatrix
   sefsef. Sefsef then gets pushed to GA for subsequent eigenvector analysis.
*/

//Entry to computing the matrix elements
int PsociHamiltonian::generateSOCIblocks(JOUTFG & iconf, JOUTFG & jconf, vector<double> & sefsef )
{
  int have_processed = 0; //guilty until proven innocent- indicates no sefsef contributions computed
  
  vector<double> vnt;
  sefsef.clear(); // Rebuild data structure from scratch ( nsefi X nsef j) in size
  
#ifdef DETAILEDHAMCHECK
  if ( iconf.completed_state != 1|| jconf.completed_state != 1) {
    cerr << "generateSOCIblocks: One of both input confs are not in completed_state: Aborting " << endl;
    GA::Terminate();
  }
#endif
  
  int ndeti = iconf.ndeti;
  int ndetj = jconf.ndeti;

  int nsefi = iconf.nsefi;
  int nsefj = jconf.nsefi;

  int sefsefSize = nsefi * nsefj;
  sefsef.resize( sefsefSize ) ;

  vector<double> detdet(ndeti * ndetj, 0);
  
//#ifdef DETAILEDHAMCHECK
//now this check is broken since the size <> actual number of valid entries
  if ( iconf.alpha_spins.size() != ndeti || jconf.alpha_spins.size() != ndetj ) {
    cerr << "Expected mismatch: ndeti not alpha_spins.size(): Aborting " << endl;
    GA::Terminate();
  }
  if ( iconf.beta_spins.size() != ndeti || jconf.beta_spins.size() != ndetj ) {
    cerr << "Expected mismatch: ndeti not beta_spins.size(): Aborting " << endl;
    GA::Terminate();
  }
//#endif
  
  //If greater than 4, no non-zero interactions

  int deltaOpn =  abs( iconf.num_singly - jconf.num_singly );
  if ( deltaOpn > 4 ) {
    return( have_processed );
  }
  
  // Generate excitation specifics
  
  vector<int> iex, jex; //Needed by single exc(nmarkb=1) + double exc only empty otherwise
  
  // Note alphas and betas are already sorted 
  
  int indexi = iconf.index;
  int indexj = jconf.index;
  
  // Generate the set of excited orbitals between these two configurations
  
  vector<int> iocc( l_nbf, 0 );
  int sumelec = 0;
  
  for( int ia=0; ia < iconf.num_doubly; ++ia ) {
    sumelec += 2;
    iocc[ iconf.doubly[ia] ] = 2;
  }
  for( int ia=0; ia < iconf.num_singly; ++ia ) {
    ++sumelec;
    iocc[ iconf.singly[ia] ] = 1;
  }

  if ( sumelec != l_nelec ) {
    cerr << "Probable corrupted ioutfg object. Sum iconf electrons != nelec " << sumelec << endl;
    GA::Terminate();
  }
  
  sumelec = 0;
  vector<int> jocc( l_nbf, 0 );

    for( int ja=0; ja < jconf.num_doubly; ++ja ) {
      sumelec += 2;
      jocc[ jconf.doubly[ja] ] = 2;
    }
    for( int ja=0; ja < jconf.num_singly; ++ja ) {
      ++sumelec;
      jocc[ jconf.singly[ja] ] = 1;
    }

    if ( sumelec != l_nelec ) {
      cerr << "Probable corrupted joutfg object. Sum jconf electrons != nelec " << sumelec << endl;
      GA::Terminate();
    }
    int jmarr;
    int sumiex=0;
    int sumjex=0;
    
//Move back to the offdiag blocks 

    // Excitation Orbital lists
    for(int i=0; i< l_nbf; ++i) {
      jmarr = iocc[i] - jocc[i];
      if ( jmarr < 0 ) {
	if ( (sumjex-jmarr) > 2 ) return( have_processed );
	for(int j=-1; j>=jmarr; --j) {
	  jex.push_back( i );
	  ++sumjex;
	}
      } else if ( jmarr > 0 ) {
	if ( (sumiex+jmarr) > 2 ) return( have_processed );
	
	for(int j=1; j<=jmarr; ++j) {
	  iex.push_back( i );
	  ++sumiex;
	}
      }
    }

#if 0
    cout << "Do we have IEX orbs ? " << endl;
    jmarr=0;
    for(int i=0; i< iex.size(); ++i) {
      cout << " IEX(+1) " << iex[i]+1 ;
    }
    cout << endl;
    for(int i=0; i< jex.size(); ++i) {
      cout << " JEX(+1) " << jex[i]+1 ;
    }
    cout << endl;
#endif

// Begin the work 

  pair<int,double> info; // For future considerations

  int status = 0;
  int sefstatus = 0;
  
  // Diagonal 
  
  if ( indexi == indexj ) {
    sefstatus = 0;
    int status = generateDiagonalBlocks( iconf, iocc, vnt, info );
    if ( status == 1 ) {
      sefstatus = generateSEFCIG(iconf, jconf, iex, jex, vnt, detdet );
    }
  }
  
  // Offdiagonal

  if ( indexi != indexj ) {
    sefstatus = 1;
    int status =  generateOffdiagBlocks( iconf, jconf, iocc, jocc, iex, jex, vnt, detdet, info ) ;
    
    if ( status == 1 ) {
      sefstatus = generateSEFCIG(iconf, jconf, iex, jex, vnt, detdet );
    }
    if ( status == 0 ) sefstatus = 0;
  }
  if ( status == 2 ||  sefstatus == 1 )  ciout(iconf, jconf, detdet, sefsef );
  
  return( sefstatus );
}
//If return <> 0 then process through sefcig.
//If return < 0 we have an untrapped error

/* This is a modified version of the Hamiltonian useful in constructing the 1-e CI density
   The method is based on the work of Seth, Wagner, Ermler, Tilson
   The implementation is loosely based on the code by M. Seth and W. Ermler named DENSITY.F
   Approximate date of May 1999.

   We are going to build a new Hamiltonian matrix ( in parallel) for use of the CIDEN methods
   This is a little overkill since the matrix has many fewer non-zero elements ( it is a 1-e matrix)
*/

// UNTESTED
int PsociHamiltonian::generateSOCIblocksForCIDEN(JOUTFG & iconf, JOUTFG & jconf, vector<double> & sefsef )
{
/*
     -- from DENSITY.f
     In the 1-particle density matrix defined in terms of orthogonal MO
     elements over different CSFs can only
     contribute to one element of the density matrix
     Corresponding to the basis function where they differ.
     This element is always 1.0 so this is done in sefcigCIDEN
*/

  int have_processed = 0; //guilty until proven innocent- indicates no sefsef contributions computed
  vector<double> vnt;
  sefsef.clear(); // Rebuild data structure from scratch ( nsefi X nsef j) in size
  
#ifdef DETAILEDHAMCHECK
  if ( iconf.completed_state != 1 || jconf.completed_state != 1) {
    GA::error("generateSOCIblocksForCIDEN: One of both input confs are not in completed_state: Aborting ",-1);
  }
#endif
  
  int ndeti = iconf.ndeti;
  int ndetj = jconf.ndeti;

  int nsefi = iconf.nsefi;
  int nsefj = jconf.nsefi;

  int sefsefSize = nsefi * nsefj;
  sefsef.clear();
  sefsef.resize( sefsefSize, 0.0 ); // ensures all are zero 

  vector<double> detdet(ndeti * ndetj, 0);
  
  //If greater than 2, no non-zero interactions

  int deltaOpn =  abs( iconf.num_singly - jconf.num_singly );
  if ( deltaOpn > 2 ) {
    return( have_processed );
  }

  if ( jconf.jsym != iconf.jsym ) {
    return( have_processed );
  } 
  
  vector<int> iex, jex; //Needed by single exc(nmarkb=1) + double exc only empty otherwise
  
  // Note alphas and betas are already sorted 
  
  int indexi = iconf.index;
  int indexj = jconf.index;
  
  if ( indexi != indexj ) {

  // Generate the set of excited orbitals between these two configurations
  
  vector<int> iocc( l_nbf, 0 );
  int sumelec = 0;
  
  for( int ia=0; ia < iconf.num_doubly; ++ia ) {
    sumelec += 2;
    iocc[ iconf.doubly[ia] ] = 2;
  }
  for( int ia=0; ia < iconf.num_singly; ++ia ) {
    ++sumelec;
    iocc[ iconf.singly[ia] ] = 1;
  }

  if ( sumelec != l_nelec ) {
    cerr << "Probable corrupted ioutfg object. Sum iconf electrons != nelec " << sumelec << endl;
    GA::Terminate();
  }
  
  sumelec = 0;
  vector<int> jocc( l_nbf, 0 );

    for( int ja=0; ja < jconf.num_doubly; ++ja ) {
      sumelec += 2;
      jocc[ jconf.doubly[ja] ] = 2;
    }
    for( int ja=0; ja < jconf.num_singly; ++ja ) {
      ++sumelec;
      jocc[ jconf.singly[ja] ] = 1;
    }

    if ( sumelec != l_nelec ) {
      cerr << "Probable corrupted joutfg object. Sum jconf electrons != nelec " << sumelec << endl;
      GA::Terminate();
    }
    int jmarr;
    int sumiex=0;
    int sumjex=0;
    
//Move back to the offdiag blocks 

    // Excitation Orbital lists
    for(int i=0; i< l_nbf; ++i) {
      jmarr = iocc[i] - jocc[i];
      if ( jmarr < 0 ) {
	if ( (sumjex-jmarr) > 1 ) return( have_processed );
	for(int j=-1; j>=jmarr; --j) {
	  jex.push_back( i );
	  ++sumjex;
	}
      } else if ( jmarr > 0 ) {
	if ( (sumiex+jmarr) > 1 ) return( have_processed );
	
	for(int j=1; j<=jmarr; ++j) {
	  iex.push_back( i );
	  ++sumiex;
	}
      }
    }

  } // if ( indexi != indexj ) 

  pair<int,double> info; // For future considerations

  int sefstatus = 0;
  
  sefstatus = generateSEFCIGCIDEN(iconf, jconf, iex, jex, detdet );
  if ( sefstatus == 1 ) ciout(iconf, jconf, detdet, sefsef );
  
/* We could transform the blocks here and now without ever going to the GA space again.
   Then the accumulated density matrix could be simply GOP'd at the end.

   It is easier to do now since we already have the iex/jex. But we now need the CI vectors
   and the number of basis functions. We can transform back to AOs in a later step.

   We already have NBF as l_nbf. Perform local block transformations and assumulate at the end
   densityMatrix[numRoots][ nbf*nbf ];

*/
  return( sefstatus );
}

// Beginnings of the SOCIG code: Previously called socig

int PsociHamiltonian::socig(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iex, vector<int> & jex, vector<double> & vnt, 
			    vector<double> & detdet )
{ 
 /*
  evaluate the spin-orbit matrix elements for configurations iconf
  and jconf (which differ by a single excitation and are not in
  maximum alignment) by evaluating the hamiltonian matrix over the
  corresponding determinants.
*/
  
  int nmarka, nmarkb;
  
  // For the given iconf/jconf pair compute nmaeka and nmarkb
  
  int ndeti = iconf.ndeti;
  int ndetj = jconf.ndeti;
  
  // Fill ifig(a/b)[nbf] with 1s to denote open shell orbs, 0s otherwise.
  // 1s are inserted in the absolute position specified by alpha(beta)_spins values.
  
#ifdef DETAILEDHAMCHECK
/*
  if ( iconf.alpha_spins.size() != iconf.beta_spins.size() ) {
    cerr << " Erroneous spin lists. alpha(ndeti) <> beta(ndeti) for index " << iconf.index;
    GA::Terminate();
  }
*/
#endif
  
  int nalphai, nbetai;
  int nalphaj, nbetaj;
  
  int jsymi = iconf.jsym;
  int jsymj = jconf.jsym;
  
  int ncount;

  for(int ii=0; ii< ndeti; ++ii ) {
    
    vector<int> aorbsi = iconf.alpha_spins[ii];
    vector<int> borbsi = iconf.beta_spins[ii];
    
    nalphai = iconf.nalpha[ii]; //aorbsi.size();
    nbetai = iconf.nbeta[ii];  //borbsi.size();
    
    vector<int> ifiga( l_nbf, 0);
    vector<int> ifigb( l_nbf, 0);
    
    for(int ka=0; ka< nalphai; ++ka) {
      ifiga[ aorbsi[ka] ] = 1;
    }
    
    for(int ka=0; ka< nbetai; ++ka) {
      ifigb[ borbsi[ka] ] = 1;
    }
    
    for( int jj =0; jj< ndetj; ++jj ) {  
      vector<int> aorbsj = jconf.alpha_spins[jj];
      vector<int> borbsj = jconf.beta_spins[jj];
      
      nalphaj = jconf.nalpha[jj]; //aorbsj.size();
      nbetaj = jconf.nbeta[jj];  //borbsj.size();
      
      nmarka=0;
      for (int ka=0; ka < nalphaj; ++ka ) {
	
	if ( ifiga[ aorbsj[ka] ] != 1 )  ++nmarka ;
	if ( nmarka > 1 ) break; 
      }
      nmarkb=nmarka;
      for (int ka=0; ka < nbetaj; ++ka ) {
	if ( ifigb[ borbsj[ka] ] != 1 ) ++nmarkb;
	if ( nmarkb > 1 ) break;
      }

      if ( nmarka <= 1 && nmarkb <= 1 ) {
/*
  single excitation
  the excited orbitals were determined in the calling routine
  
  iex(1)       the excited orbital of configuration i
  jex(1)       the excited orbital of configuration j
  
  determine the sign (ncount is the number of sign changes)
*/
	ncount = iconf.phases[ii] + jconf.phases[jj];

	if ( nalphai < nalphaj ) {
/*
   i-beta j-alpha case
*/
	  for( int i=0; i< nbetai; ++i) {
	    if (  borbsi[i] < iex[0] ) ++ncount;
	  }
	  for( int i=0; i< nalphaj; ++i) {
	    if (  aorbsj[i] > jex[0] ) ++ncount;
	  }
	  
// ly or lx ?
	  if ( (jsymi+jsymj)%2 == 0 ) { // ly elements
	    ++ncount;
	  } else if( (jsymi%2) == 1 ) { // lx elements
	    ++ncount;
	  }
	} else if ( nalphai > nalphaj ) {
/*
       # i-alpha j-beta case
*/
	  for( int i=0; i< nalphai; ++i) {
	    if (  aorbsi[i] > iex[0] ) ++ncount;
	  }
	  for( int i=0; i< nbetaj; ++i) {
	    if (  borbsj[i] < jex[0] ) ++ncount;
	  }
	  
// ly or lx ? 
	  if ( (jsymi+jsymj)%2 != 0 ) { // lx elements
	    if ( (jsymi%2) == 1 ) ++ncount;            
	  }
	} else {
/*
       # i-alpha j-alpha and i-beta j-beta cases
*/
	  if ( nmarka == 0 ) ++ncount;
	  if ( (jsymi%2) == 1 ) ++ncount;
	  
	  if( abs(iex[0]-jex[0]) > 1 ) {
	    int kbs = max( iex[0],jex[0] ) - 1;
	    int ka = iex[0] + jex[0] - kbs;
	    if ( nmarka != 0 ) {
	      for(int ia=0; ia< nalphai; ++ia) {
		if ( aorbsi[ia] > kbs ) break;
		if ( aorbsi[ia] >= ka ) ++ncount;
	      }
	    } else {
	      for(int ia=0; ia< nbetai; ++ia) {
		if ( borbsi[ia] > kbs ) break;
	        if ( borbsi[ia] >= ka ) ++ncount;
	      }
	    }
	  } // End conditions
	} // nalphai nalphaj comparison
	int index = ii*ndetj + jj;
	( ncount % 2 == 1 ) ? detdet[index] = -vnt[0] : detdet[index] = vnt[0];
      } // nmarka nmarkb comparison
    } //jj
  } //ii
  
  return( 0 );
}

int PsociHamiltonian::generateSEFCIG(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iex, vector<int> & jex, vector<double> & vnt, vector<double> & detdet )
{
/*
  evaluate the kinetic-energy, coulomb-interaction matrix-element
  block over configurations iconf and jconf (which differ by a single
  or double excitation and are not in maximum alignment) by
  evaluating matrix elements over the corresponding determinants.
  the determinants must have the same spin projection to have
  non-zero matrix elements.
*/

// DO NOT empty detdet: socig may have populated some of it

  int status = 0; //guilty until proven innocent
  
  int ndeti = iconf.ndeti;
  int ndetj = jconf.ndeti;
  
#ifdef DETAILEDHAMCHECK
  if (detdet.size() != ndeti*ndetj ) {
    cerr << "detdet has wrong length " << detdet.size() << " ndeti*ndetj is " << ndeti*ndetj << endl;
    GA::Terminate();
  }
#endif

  int msi, msj; // spin projections
  
  int nalphai, nbetai;
  int nalphaj, nbetaj;
  int nopeni;
  
  int index; // for detdet index
  
#if 0
  cout << "SEFCIG i and j are " << iconf.index << " " << jconf.index << endl;
  cout << "ndeti is " << ndeti << " ndetj is " << ndetj <<  endl;
  cout << "TEST a/borbs " << endl;
//    for(int i=0; i< ndeti; ++i) {
//    cout << " " << iconf.alpha_spins[i] << " " <<  iconf.beta_spins[i] << endl;
//  }
#endif

  for(int ii=0; ii< ndeti; ++ii ) {
    
    vector<int> aorbsi = iconf.alpha_spins[ii];
    vector<int> borbsi = iconf.beta_spins[ii];
    
    nalphai = iconf.nalpha[ii]; //aorbsi.size();
    nbetai = iconf.nbeta[ii];  //borbsi.size();
    
    vector<int> ifiga( l_nbf, 0);
    vector<int> ifigb( l_nbf, 0);
    
    for(int ka=0; ka< nalphai; ++ka) { //iconf-centric data structure
      ifiga[ aorbsi[ka] ] = 1; 
    }
    for(int ka=0; ka< nbetai; ++ka) { // ditto
      ifigb[ borbsi[ka] ] = 1;
    }
    
    nopeni = iconf.num_singly;
    msi = iconf.spin_projection[ii];
    
    int nmarka, nmarkb;
    int nij[4];
    
    for( int jj =0; jj< ndetj; ++jj ) {
      index = ii*ndetj + jj;
      msj = jconf.spin_projection[jj];
      if ( msi == msj ) {  
	vector<int> aorbsj = jconf.alpha_spins[jj];
	vector<int> borbsj = jconf.beta_spins[jj];
	
         nalphaj = jconf.nalpha[jj]; //aorbsi.size();
         nbetaj = jconf.nbeta[jj];  //borbsi.size();
	
	//nopenj = jconf.num_singly;
	nmarka=0;
	for (int ka=0; ka < nalphaj; ++ka ) {
	  //cout << ii<<" "<<jj<<" "<<"ALPHA " << endl;
	  if ( ifiga[ aorbsj[ka] ] != 1 ) {
            //cout << "a nonmatched alpha is " << aorbsj[ka] << endl;
	    ++nmarka;
	    if ( nmarka > 2 ) {
	      break;
	    }
	    nij[nmarka-1] = ka;
          }
	}
	nmarkb=nmarka;
	for (int ka=0; ka < nbetaj; ++ka ) {
	  //cout << ii<<" "<<jj<<" ""BETA " << endl;
	  if ( ifigb[ borbsj[ka] ] != 1 ) {
            //cout << "a nonmatched beta is " << borbsj[ka] << endl;
	    ++nmarkb;
	    if ( nmarkb > 2 ) { 
	      break;
	    }
	    nij[nmarkb-1] = ka;
          }
	}
	
	//cout << "nmarks are " << nmarka<<" "<<nmarkb<<endl;
	if ( nmarkb <=2 && nmarka <=2 ) {
/*
     Only non trivial scenarios are:
     nmarka, nmarkb combinations:
                    (0,0); (0,1), (1,1); (0,2), (1,2), (2,2)
*/
// No excitation
	  
	  if (nmarkb == 0 ) {
	    int iorb, jorb;
	    int mm = 0;
	    detdet[index] = vnt[mm];
	    
	    if ( nopeni > 1 ) {
	      for(int iopen=1; iopen < nopeni; ++iopen) {
		iorb = iconf.singly[iopen];
		for(int jopen=0; jopen < iopen; ++jopen) {
		  jorb = jconf.singly[jopen];
		  ++mm;
		  if( ifiga[iorb] == ifiga[jorb] ) detdet[index] -= vnt[mm];
		}
	      }
	    }
	    status = 1;
	  } else if( nmarkb == 1 ) {
/*
  single excitation

  the excited orbitals were determined in the calling routine

     iex(1)       the excited orbital of configuration i
     jex(1)       the excited orbital of configuration j
*/

	    int ncount = iconf.phases[ii] + jconf.phases[jj];
	    if ( abs( iex[0] - jex[0] ) > 1 ) {
	      int kbs = max( iex[0],jex[0] )-1;
	      int ka = iex[0]+jex[0]-kbs;
	      
	      if ( nmarka != 0 ) {
		for( int ia=0; ia < nalphai; ++ia) {
		  if( aorbsi[ ia ] > kbs ) break;
		  if( aorbsi[ ia ] >= ka ) ++ncount;
		}
	      } else {
		for( int ia=0; ia < nbetai; ++ia) {
		  if( borbsi[ ia ] > kbs ) break;
		  if( borbsi[ ia ] >= ka ) ++ncount;
		}
	      }
	    }
	    
	    detdet[index] = vnt[0]; 
	    
	    int temp;
	    for(int ia=0; ia < l_nbf; ++ia ) {
	      ( nmarka == 0 )? temp = ifigb[ ia ]: temp =  ifiga[ ia ];
	      
	      if ( temp != 0 ) {
		if ( ia != iex[0] && ia != jex[0] && ifiga[ia] + ifigb[ia] != 2 ) { 
		  detdet[index] -= l_integrals->fetch2eIntegralValue(iex[0],ia,ia,jex[0] );
		}  
	      }
	    } // end ia
	    if( ncount%2 == 1 ) detdet[index] = -detdet[index];
	    
	    status = 1;
	  } else if ( nmarka == 1 ) {
 /*
  double excitation (product of single excitations from the
  alpha and beta sets)

  determine the excited orbitals
     iex(2) and iex(3)      the excited orbitals of the alpha set
     jex(2) and jex(3)      the excited orbitals of the beta set
*/

	    int iex2, jex2, iex3, jex3; 
	    int ncount = iconf.phases[ii] + jconf.phases[jj];
	    int nij1 = nij[0];
	    int nij2 = nij[1];
	    iex2 = aorbsj[ nij1 ];
	    jex2 = borbsj[ nij2 ]; 
	    
	    iex3 = iconf.sumsa[ ii ] - jconf.sumsa[ jj ] + iex2;
	    jex3 = iconf.sumsb[ ii ] - jconf.sumsb[ jj ] + jex2;
	    
	    //cerr << "iconf index " << iconf.index<<" " << jconf.index<<endl;
	    //cerr <<"Entry to again "<< iex2<<" "<<jex2<<" "<<iex3<<" "<<jex3<<endl;
	    //cerr <<"SUMS"<<" "<<iconf.sumsa[ ii ] <<" "<<iconf.sumsb[ ii ]<<" "<<jconf.sumsa[ jj ]<<" "<<jconf.sumsb[ jj ]<<endl;
	    
/*
  determine the signs (ncount is the number of sign changes)
*/

	    if ( abs( iex2 - iex3) > 1 ) { 
	      int kbs = max(iex2,iex3)-1;
	      int ka = iex2+iex3-kbs;
	      for(int ia=0; ia< nalphai; ++ia) {
		if( aorbsi[ ia ] > kbs ) break;
		if( aorbsi[ ia ] >= ka ) ++ncount;
	      }
	    }
	    if ( abs( jex2 - jex3) > 1 ) {
	      int kbs = max(jex2,jex3)-1;
	      int ka = jex2+jex3-kbs;
	      for(int ia=0; ia< nbetai; ++ia) {
		if( borbsi[ ia ] > kbs ) break;
		if( borbsi[ ia ] >= ka ) ++ncount;
	      }
	    }
	    detdet[index] = l_integrals->fetch2eIntegralValue(jex2,jex3,iex2,iex3 );
	    if( ncount%2  == 1 ) detdet[index] = -detdet[index];
            
	    status = 1;
	    
	  } else {
/*
  double excitation
 
  the excited orbitals were determined in the calling routine

     iex(1) and iex(2)      the excited orbitals of configuration i
     jex(1) and jex(2)      the excited orbitals of configuration j

  determine the sign (ncount is the number of sign changes)
*/
	    int ncount = iconf.phases[ii] + jconf.phases[jj] + nij[1] - nij[0] -1;
	    
	    if ( iex[1] - iex[0] > 1 ) {
	      if ( nmarka != 0 ) {
		for(int ia=0; ia< nalphai; ++ia) {
		  if( aorbsi[ ia ] >= iex[1] ) break;
		  if( aorbsi[ ia ] > iex[0] ) ++ncount;
		}
	      } else {
		for(int ia=0; ia< nbetai; ++ia) {
		  if( borbsi[ ia ] >= iex[1] ) break;
		  if( borbsi[ ia ] > iex[0] ) ++ncount;
		}
	      }
	    }
	    detdet[index] = vnt[0] - vnt[1];
	    if( ncount%2  == 1 ) detdet[index] = -detdet[index];
	    
	    status = 1;
	    
	  } // nmarka and nmarkb <=2
	}
      }// msi <> msj
    }// next jj 
  }// next ii 
  
  status = 1; // Not too useful, eh
  return( status );
}
/* A pruned version of SEFCIG that only captures interactions relevant to CIDEN

   Experimental - not functional
*/
int PsociHamiltonian::generateSEFCIGCIDEN(JOUTFG & iconf, JOUTFG & jconf, vector<int> & iex, vector<int> & jex, vector<double> & detdet )
{
/*
  evaluate the kinetic-energy, coulomb-interaction matrix-element
  block over configurations iconf and jconf (which differ by a single
  or double excitation and are not in maximum alignment) by
  evaluating matrix elements over the corresponding determinants.
  the determinants must have the same spin projection to have
  non-zero matrix elements.
*/

// DO NOT empty detdet: socig may have populated some of it

  int status = 0; //guilty until proven innocent
  
  int ndeti = iconf.ndeti;
  int ndetj = jconf.ndeti;
  
#ifdef DETAILEDHAMCHECK
  if (detdet.size() != ndeti*ndetj ) {
    cerr << "CIDEN detdet has wrong length " << detdet.size() << " ndeti*ndetj is " << ndeti*ndetj << endl;
    GA::Terminate();
  }
#endif

  int msi, msj; // spin projections

  double dZero = 0.0;
  double dOne = 1.0;
  
  int nalphai, nbetai;
  int nalphaj, nbetaj;
  
  int index; // for detdet index
  
#if 0
  cout << "SEFCIG i and j are " << iconf.index << " " << jconf.index << endl;
  cout << "ndeti is " << ndeti << " ndetj is " << ndetj <<  endl;
  cout << "TEST a/borbs " << endl;
//    for(int i=0; i< ndeti; ++i) {
//    cout << " " << iconf.alpha_spins[i] << " " <<  iconf.beta_spins[i] << endl;
//  }
#endif

  for(int ii=0; ii< ndeti; ++ii ) {
    
    vector<int> aorbsi = iconf.alpha_spins[ii];
    vector<int> borbsi = iconf.beta_spins[ii];
    
    nalphai = iconf.nalpha[ii]; //aorbsi.size();
    nbetai = iconf.nbeta[ii];  //borbsi.size();
    
    vector<int> ifiga( l_nbf, 0);
    vector<int> ifigb( l_nbf, 0);
    
    for(int ka=0; ka< nalphai; ++ka) { //iconf-centric data structure
      ifiga[ aorbsi[ka] ] = 1; 
    }
    for(int ka=0; ka< nbetai; ++ka) { // ditto
      ifigb[ borbsi[ka] ] = 1;
    }
    
    //nopeni = iconf.num_singly;
    msi = iconf.spin_projection[ii];
    
    int nmarka, nmarkb;
    //int nij[4]; //probably can make this a TWO
    
    for( int jj =0; jj< ndetj; ++jj ) {
      index = ii*ndetj + jj;
      msj = jconf.spin_projection[jj];
      if ( msi == msj ) {  
	vector<int> aorbsj = jconf.alpha_spins[jj];
	vector<int> borbsj = jconf.beta_spins[jj];
	
         nalphaj = jconf.nalpha[jj]; //aorbsi.size();
         nbetaj = jconf.nbeta[jj];  //borbsi.size();
	
	//nopenj = jconf.num_singly;
	nmarka=0;
	for (int ka=0; ka < nalphaj; ++ka ) {
	  //cout << ii<<" "<<jj<<" "<<"ALPHA " << endl;
	  if ( ifiga[ aorbsj[ka] ] != 1 ) {
            //cout << "a nonmatched alpha is " << aorbsj[ka] << endl;
	    ++nmarka;
	    if ( nmarka > 1 ) {
	      break;
	    }
	    //nij[nmarka-1] = ka;
          }
	}
	nmarkb=nmarka;
	for (int ka=0; ka < nbetaj; ++ka ) {
	  //cout << ii<<" "<<jj<<" ""BETA " << endl;
	  if ( ifigb[ borbsj[ka] ] != 1 ) {
            //cout << "a nonmatched beta is " << borbsj[ka] << endl;
	    ++nmarkb;
	    if ( nmarkb > 1 ) { 
	      break;
	    }
	    //nij[nmarkb-1] = ka;
          }
	}
	
	//cout << "nmarks are " << nmarka<<" "<<nmarkb<<endl;
	
        if ( nmarka > 1 || nmarkb > 1 ) {
	  detdet[index] = dZero;
	  
        } else if ( nmarkb == 0 ) {
          //cout << "at the nmarkb=0"<<endl;
	  detdet[index] = dOne;
	  
        } else {
	  
          //cout << "at the the nmarb0 or nmarka0,1"<<endl;
	  int ncount = iconf.phases[ii] + jconf.phases[jj];
	  if ( abs( iex[0] - jex[0] ) > 1 ) {
	    int kbs = max( iex[0],jex[0] )-1;
	    int ka = iex[0]+jex[0]-kbs;
	    
	    if ( nmarka != 0 ) {
	      for( int ia=0; ia < nalphai; ++ia) {
		if( aorbsi[ ia ] > kbs ) break;
		if( aorbsi[ ia ] >= ka ) ++ncount;
	      }
	    } else {
	      for( int ia=0; ia < nbetai; ++ia) {
		if( borbsi[ ia ] > kbs ) break;
		if( borbsi[ ia ] >= ka ) ++ncount;
	      }
	    }
	  }
	  detdet[index] = dOne;
	  if( ncount%2 == 1 ) detdet[index] = -detdet[index];
        }
      }// msi <> msj
       
      // works Okay with a small test problem
      // cout << "CIDEN ii jj index detdet "<<ii<<" "<<jj<<" "<<index<<" "<<detdet[index]<<endl;

    }// next jj 
  }// next ii 
  
  status = 1; // Not too useful, eh
  return( status );
}

// Local call: Mostly for debugging purposes
void PsociHamiltonian::print2Darray(int nrow, int ncol, vector<double> & detdet )
{
  cout << "A printout of block requested " << endl;
  cout << "Rows are " << nrow << " Cols are " << ncol << endl;
  cout << endl;
  
  if ( detdet.size() != nrow * ncol ) {
    cerr << "Error in print2Darray: detdet and i-j mismatch " << endl;
    GA::Terminate();
  }
  
  int irow = 0;
  for(int i=0; i< nrow; ++i ) {
    cout << "next row = " << ++irow << " " << endl;
    for( int j=0; j< ncol; ++j ) {
     cout << std::setprecision (15) << detdet[ (ncol*i)+j ] << " ";
    } 
    cout << endl;
  }
  return;
}

// This is purely determinant related data but only needed by hamiltonian construction. So put it here.
/* This is a significant change from the original. Here we simply grab a linked list of 
   pair of dets/sefs that couple. We figure out phase (sign) and pairs in subsequent transforms

   A(pair) == A(maxsef, 0:1); where A(,0) first det index and (,1) is the second
*/

// maxdet and maxsef are maximums PER SPATIAL and not global sums
// Treat the arrays as having dimensions of A( maxsef, maxdet );

int PsociHamiltonian::determinantToDoublegroupMaps( vector<pair<int,int> > & coef ) 
{
  if ( l_maxdet == 0 || l_maxsef == 0 ) {
    cerr << "determinantToDoublegroupMaps: maxdet or maxsef not set " << endl;
    GA::Terminate();
  }
  
  coef.clear();
  
  int ij = 0;
  
  //  if ( l_nelec % 2 == 0 ) { go ahead and create regardless of nelec 
  for(int ii=0; ii < l_maxsef; ++ii) {
    coef.push_back(  pair<int,int>( ij, ij+1) ) ;
    ij += (l_maxdet + 2 );
  }
  //  }
  return(0) ;
}

//Timer wrapped version
void  PsociHamiltonian::transformCDCt( int ncpi, int ncpj, int nopeni, int nopenj, int nsefi, int ndeti, int nsefj, int ndetj, int & rows, int & cols, vector<double>  & sefsef, 
				       vector<pair<int,int> > & coef, vector<double> & detdet, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  transformCDCt( ncpi, ncpj, nopeni, nopenj, nsefi, ndeti, nsefj, ndetj, rows, cols, sefsef, coef, detdet);
  info.second += psociTime() - timein; //Accumulate times
  info.first = GA::nodeid();
  return;
}

/*
    Basic dgemm wrapper: m,n,k are likely to mostly be 8 to 16 in value increasing to as many as
    1024.

    C( m, n ) = C( m, n ) + A( m, k ) * B( n, k )

    A arrives as a packed array: A( m, pair<int,int> ); pair reveals the 2 non-zero terms
    
    ncp ==1 summation, ncp ==2 difference
    No support for transposition is  provided
*/

// Local call

void  PsociHamiltonian::transformCDCt( int ncpi, int ncpj, int nopeni, int nopenj, int nsefi, int ndeti, int nsefj, int ndetj, int & rows, int & cols, vector<double>  & sefsef, 
				       vector<pair<int,int> > & coef, vector<double> & detdet )
{
#ifdef DETAILEDHAMCHECK
  if( coef.size() < nsefj ) {
    cerr << " transformCDCt: coef of unexpected length: size is " << coef.size() << endl;
    GA::Terminate();
  }
  if ( detdet.size() != ndetj * ndeti ) {
    cerr << " transformCDCt: detdet of unexpected size: size is " << detdet.size() << endl;
    GA::Terminate();
  }
#endif
  
#ifdef BLAS
  cout << " Should probably use a good dgemm. But for now keep it simple. " << endl;
  cout << " Use naive approach for now " << endl;
#endif

  double rthlf = sqrt ( 0.5 );
  
  double temp;
  
  // Perform first haftTransform: C = C + dtsefj * detdet
  // detdet(ndeti,ndetj) * coef (nsefj, pair<int,int>) -> detsef( ndeti, nsefj )
  
  int ka,kb;
  
  double multi = 1.0; 
  if ( ncpj == 2 ) multi =-1.0;
  
  if ( l_nelec % 2 == 1 ) {
    sefsef = detdet;
     rows = ndeti;
     cols = ndetj;
     return;
  }
  
  int indDetsef, indDetdet; 
  int offset;
  
  vector<double> detsef;
  
  if ( nopenj != 0 ) { 
    detsef.resize( ndeti * nsefj);
    
    for(int ii=0; ii < ndeti; ++ii ) {
      indDetsef = ii * nsefj;
      indDetdet = ii * ndetj;
      for(int jj=0; jj < nsefj; ++jj ) {
	ka = coef[jj].first;
	kb = coef[jj].second;
	offset = jj * l_maxdet;
	temp = detdet[ indDetdet + (ka - offset) ];
	temp += multi * detdet[indDetdet + (kb - offset) ];
	detsef[ indDetsef + jj ] = temp * rthlf;
      }
    }
  } else {
    detsef = detdet;
  }

#if 0
  cout << " nopenj is " << nopenj << " Output a DETDET array " << endl;
  print2Darray( ndeti, ndetj, detdet );

  cout << " nopenj is " << nopenj << " Output a DETSEF array " << endl;
  print2Darray( ndeti, nsefj, detsef );
#endif

  // Perform SecondHalfTransform: sefsef(nsefi, nsefj) = coef(nsefi, ndeti) * detsef( ndeti, nsefj )
  // detsef(ndeti,nsefj) * coef (nsefi, pair<int,int>) -> sefsef( nseti, nsefj )
  
  multi = 1.0;
  if ( ncpi == 2 ) multi =-1.0;
  
  int indSefsef;
  if ( nopeni != 0 ) {
    sefsef.resize( nsefi * nsefj );
    rows = nsefi;
    cols = nsefj;
    
    for(int ii=0; ii < nsefi; ++ii ) {
      indSefsef = ii * nsefj; 
      ka = coef[ii].first;
      kb = coef[ii].second; // These appear correct
      offset = ii * l_maxdet; // coef numbers are based on the fullsize matrix
      for(int jj=0; jj < nsefj; ++jj ) {
	temp = rthlf * detsef[  (ka - offset )  * nsefj + jj ];
	temp += multi * rthlf * detsef[ (kb - offset ) * nsefj + jj ];
	sefsef[ indSefsef + jj ] = temp;
      }
    }
  } else {
    sefsef = detsef;
    rows = ndeti;
    cols = nsefj;
  }
#if 0
  cout << "Output a final SEFSEF array for  " << rows<<" "<<cols<<endl;
  print2Darray( rows, cols, sefsef );
#endif
  
  return;
}

//
void PsociHamiltonian::ciout( JOUTFG & iconf, JOUTFG & jconf, vector<double> & detdet, vector<double> & sefsef )
{
  /* Transform detdet matrix double group subblock. Then determine actual GA Hamiltonian
     columns for insertion
  */
  int ndeti = iconf.ndeti;
  int nsefi = iconf.nsefi;
  
  int ndetj = jconf.ndeti;
  int nsefj = jconf.nsefi;
  
  int isym = iconf.jsym;
  int jsym = jconf.jsym;
  
#ifdef DETAILEDHAMCHECK
  if ( l_ksym == 0 ) {
    cerr << "ciout: ksym is probably not specified: Must set with void PsociHamiltonian::setKsym " << endl;
    GA::Terminate();
  }
#endif
  
  int ksym = l_ksym; // Same for all 
  
  int nopeni = iconf.num_singly;
  int nopenj = jconf.num_singly;
  
  vector<pair<int,int> > coef;
  
  determinantToDoublegroupMaps( coef ); //This should be a static map so we do not repeatedly compute this
  
  // Choose which coefs (a+b or a-b) I and J require. Based on symmetry
  int ncpi, ncpj;
  if ( l_nelec % 2 == 1 ) {
    ncpi = 2;
    ncpj = 2;
  } else {
    ( isym % 2 != ksym % 2 ) ? ncpi = 2: ncpi = 1;
    ( jsym % 2 != ksym % 2 ) ? ncpj = 2: ncpj = 1;
  }
  
/* Perform the C*D*C*t similarity transform. I am not using dgemm codes for now
     because of the sparsity (2 values) in the C terms.
*/

// all checks for mel and openi/j are now within CDCt
// Final result is always in sefsef

  //sefsef size specified in transform
  
  int rows;
  int cols;
  transformCDCt( ncpi, ncpj, nopeni, nopenj, nsefi, ndeti, nsefj, ndetj, rows, cols, sefsef, coef, detdet);
  
/* At this point we have the full sefsef hamiltonian between configurations I and J
   Now we need to package this up for sending to GA 
 ( and eventually to DRA is need be )
*/
  
  return;
}

