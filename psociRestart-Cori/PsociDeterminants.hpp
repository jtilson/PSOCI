/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

* New implementation of PsociDeterminants:

 Classes: 

 Description: 

 History:

 July 1, 2011: I added back to the joutfg type the variables nbf,nelec,mx_openshells,ksym. Not because
               we actually need them. After fetching out of GA, these values get lost. Now they are NEVER
               subsequently used BUT it is very confusing to the user if one ir printing out these objects
               with the lost terms

**************************************************************************************/
/**
 *   @file PsociDeterminants.hpp
 *
 */

#ifndef PSOCIDETERMINANTS_H
#define PSOCIDETERMINANTS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <map>
#include <vector>

#include "ga++.h"

#include "PsociConfigs.hpp"

using namespace std;

// DANGER must make consistent with PsociGADeterminants.hpp
#ifdef USESHORT
typedef short COMPRESS_SIZE;
#else
typedef int COMPRESS_SIZE; //GA doesn't support shorts at this time
#endif

/* We could use SHORTs or smaller for many of these entries BUT we eventualy want to increase substantially the number
   of basis functions and orbs. So let's not worry about space yet. Besides, this allocation is transient. I.e., 
   it will copied into GA and deleted from this - PRIOR to hamiltonian construction
*/

/* CHANGED BOOL TO SHORT because of problem in intel composer 117 compiler
*/
struct JOUTFG {
  short completed_state;
  int index; //relative to the multimap container
  int mx_openshells; // Maximum number of openshells in the problems (formally mxocc );
  int nbf; //Number of basis ftns in the problem
  int jsym; //symmetry of spatial 
  int ksym; //State symmetry
  int nelec; // number electrons for this conf ser
  int num_doubly; // number of doubly occupied orbitals
  int num_singly; // number of singly occupied orbitals
  string config; // Grabbed from the multimap
  int nsefi; // == 4**( (nopen-1)/2 ): num sefs for configuration(spatial) i
  int ndeti; // == 2**max( nopen-1,0): num dets for configuration(spatial) i
  vector<int> singly;
  vector<int> doubly;
  vector<int> spin_projection; // spin-projection(ndeti) for the current determinant (formally named ms). = (mxopn+2*ms)/2
  vector<int> phases; // Phases for each ndeti (0,1)
  vector<int> sumsa; // Sum of alpha orbitals 
  vector<int> sumsb; // Sum of beta orbitals 
  
  vector<int> nalpha; // New entries to minimize the need to resizin on the fetchAndUnpack code
  vector<int> nbeta;  // New entries .. ditto..
  
  vector<vector<int> > alpha_spins; // alpha-spin orbs per determinant (formally called ifiga ) A[nalpha][ndeti]
  vector<vector<int> > beta_spins; //  beta-spin orbs per determinant (formally called ifigb )  B[nbeta][ndeti]
};

const int MAX_LINES_PRINT = 1000; 

const int PSOCI_ONE_OPEN=1, PSOCI_THREE_OPEN=3, PSOCI_FIVE_OPEN=5, PSOCI_SEVEN_OPEN=7;
const int PSOCI_NINE_OPEN=9, PSOCI_ELEVEN_OPEN=11;
const int PSOCI_MAX_OPEN =  PSOCI_ELEVEN_OPEN;

const int PSOCI_ZERO_OPEN=0, PSOCI_TWO_OPEN=2, PSOCI_FOUR_OPEN=4, PSOCI_SIX_OPEN=6;
const int PSOCI_EIGHT_OPEN=8, PSOCI_TEN_OPEN=10;
const int PSOCI_WHOLE_MAX_OPEN =  PSOCI_TEN_OPEN;

class PsociDeterminants {
  
private:
  PsociConfigs * local_configs;
  
  int l_mxopenshells; // Maximum number of openshells in the problems (formally mxocc );
  int l_nbf; //Number of basis ftns in the problem
  int l_ksym; //State symmetry
  int l_nelec; // number electrons 
  int l_maxspatials; // total spatials processed into determinants or not
  
  vector<JOUTFG> joutfg; //Possibly only containing a subset of processed configs
  
  bool printDeterminants;
  
  int generateDeterminants( JOUTFG & ioutfg );
  int generateSpinFractionalDeterminants( JOUTFG & ioutfg, int isgn );
  int generateEVENSpinDeterminants( JOUTFG & ioutfg );
  int generateODDSpinDeterminants( JOUTFG & ioutfg );

  int local_numcfgs; // Number ( subset) of the read spatials I am processing - could be 0;
  int global_numcfgs; // Sum across all nodes

  int local_maxdet_spatial;
  int global_maxdet_spatial;

  int local_maxsef_spatial;
  int global_maxsef_spatial;

  int local_maxdet_perspatial;
  int global_maxdet_perspatial;

  int local_maxsef_perspatial;
  int global_maxsef_perspatial;


  //vector<pair<int,string> > selcfgs; 
  
// Begin fractional electron wavefunction assembly
  int oneOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn );
  int threeOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn );
  int fiveOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn );
  int sevenOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn );
  int nineOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn );
  int elevenOpenShellFractionalDeterminants( JOUTFG & ioutfg, int isgn );
  //Begin even electron wavefunction assembly
  int zeroOpenShellEvenDeterminants( JOUTFG & ioutfg );
  int twoOpenShellEvenDeterminants( JOUTFG & ioutfg );
  int fourOpenShellEvenDeterminants( JOUTFG & ioutfg );
  int sixOpenShellEvenDeterminants( JOUTFG & ioutfg );
  int eightOpenShellEvenDeterminants( JOUTFG & ioutfg );
  int tenOpenShellEvenDeterminants( JOUTFG & ioutfg );
  //Begin odd electron wavefunction assembly
  int zeroOpenShellOddDeterminants( JOUTFG & ioutfg );
  int twoOpenShellOddDeterminants( JOUTFG & ioutfg );
  int fourOpenShellOddDeterminants( JOUTFG & ioutfg );
  int sixOpenShellOddDeterminants( JOUTFG & ioutfg );
  int eightOpenShellOddDeterminants( JOUTFG & ioutfg );
  int tenOpenShellOddDeterminants( JOUTFG & ioutfg );
  
  int computePhaseAndSort(  JOUTFG & ioutfg );
  int pushBackToJoutfg( JOUTFG & ioutfg, vector< vector<int> > & ifiga , vector< vector<int> > & ifigb, vector<int> &  ms );
  int pushBackToJoutfg( JOUTFG & ioutfg, vector< vector<int> > & ifiga , vector< vector<int> > & ifigb, vector<int> &  ms, int isign );
  
public:
  PsociDeterminants( PsociConfigs * configs );
  void printFinalDeterminants();
  void printFilename();
  void printParams();
  void printSpatials();
  int maxSpatials();
  int localSpatials();
  int numBasisFunctions();
  int fetchConfigs( int start, int stop );
  int findParams();
  int fetchConfigs( int start, int stop, pair<int,double> & info );
  int specialFetchConfigs( int start, int stop );
  int computeConfigs( int start, int stop, vector<JOUTFG> & toutfg );
  int computeConfigs( int index, JOUTFG & toutfg );
  int computeConfigs( vector<COMPRESS_SIZE> & orbmap, JOUTFG & toutfg ); // for orbmap DIRECT methods
  void printSpatialsData( int iconf );
  void printSpatialsData( JOUTFG & ioutfg );
  void tearDownDeterminants();
  int fetchLocalMaxDet();
  int fetchGlobalMaxDet();
  int fetchLocalMaxSef();
  int fetchGlobalMaxSef();
  int fetchSpinParity();
  int fetchMaxOpenShells();
  int fetchNumElectrons();
  int fetchNumBasisFunctions();
  void assembleGlobalMaxDetSef();
  void assembleGlobalNumSpatials();
  int fetchGlobalMaxSpatials();
  int fetchLocalMaxDetPerSpatial();
  int fetchGlobalMaxDetPerSpatial();
  int fetchLocalMaxSefPerSpatial();
  int fetchGlobalMaxSefPerSpatial();
  
  int fetchGlobalSymmetry();
  void assembleGlobalJobParameters();
  vector<JOUTFG> * fetchJoutfg();
  int distributeParameters();
  void copyInternalMapToVector();
  void setLocalNumcfgs( int numcfgs );
  void overrideGlobalMaxSpatials( int mxspatials );
  int specialDistributedFetchConfigs();
  int specialDistributedFetchConfigs( vector<char> & replicated_phases );
  void printLocalSpatialsData();
  
};

#endif //PSOCIDETERMINANTS_H
