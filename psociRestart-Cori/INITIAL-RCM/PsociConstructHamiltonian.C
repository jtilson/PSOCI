/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A series of functions to collect the hamiltonian restart/construction/dump
              phases of the calculation.

 History:

 

**************************************************************************************/
/**
 *   @file PsociConstructHamiltonian.C
 *
 */

#include "PsociConstructHamiltonian.hpp"

using namespace std;

void InitializeDRAservice( int g_rank, PsociDRAservices & DRAservice )
{
  DRAservice.initDRA( g_rank );
  DRAservice.reportDRAIOservers();
  DRAservice.reportDRAParams();
#if 0
  cout << "Got the DRA initialized " << g_rank << endl;
#endif
}

void InitializeDRAservice( int g_rank, PsociDRAservices & DRAservice, int IOservers )
{
  DRAservice.initDRA( g_rank, IOservers );
  DRAservice.reportDRAIOservers();
  DRAservice.reportDRAParams();
#if 0
  cout << "Got the DRA initialized " << g_rank << endl;
#endif
}

void ExitDRASERVICE( PsociDRAservices & DRAservice) 
{
  DRAservice.exitDRA();
} 

/** Allocate the four DRA data structures required to STORE GA to disk
 */
void AllocateDRA( HSOURCE state, int ndims, int maxsef, int maxsparse,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef)
{
  if ( state != CONSTRUCT ) {
    cout << "AllocateDRA: apparent mismatch: state = " << state;
  }
  
  if (GA::nodeid() == 0 ) cout << "Preparing to Allocate DRAS  Hamiltonian " << endl;
  
  pair<dra_size_t,dra_size_t> ddims;
  pair<dra_size_t,dra_size_t> rdims;
  
  rdims.first=-1;
  rdims.second=-1;
  
  d_cimat.reportName();
  d_cimat.reportFilename();
  d_cimat.setType( C_DBL );
  d_cimat.printType();
  d_cimat.setMode( DRA_RW );
  ddims.first=maxsef;
  ddims.second=maxsparse;
  d_cimat.setDimensions( ndims, ddims, rdims );
  if ( d_cimat.createDRA() != 0 )
    {
      cout << "failed to create the d_cimat: Aborting app " << endl;
      GA::error(" failed to create the d_cimat: Aborting app", GA::nodeid() );
    }
  d_cimat.printInternals();
  
  d_icicol.reportName();
  d_icicol.reportFilename();
  d_icicol.setType( C_INT );
  d_icicol.printType();
  d_icicol.setMode( DRA_RW );
  d_icicol.setDimensions( ndims, ddims, rdims );
  if ( d_icicol.createDRA() != 0 )
    {
      cout << "failed to create the d_icicol: Aborting app " << endl;
      GA::error(" failed to create the d_icicol: Aborting app", GA::nodeid() );
    }
  d_icicol.printInternals();
  
  d_number.reportName();
  d_number.reportFilename();
  d_number.setType( C_INT );
  d_number.printType();
  d_number.setMode( DRA_RW );
  ddims.first=maxsef;
  ddims.second=1;
  d_number.setDimensions( ndims, ddims, rdims );
  if ( d_number.createDRA() != 0 )
    {
      cout << "failed to create the d_number: Aborting app " << endl;
      GA::error(" failed to create the d_number: Aborting app", GA::nodeid() );
    }
  d_number.printInternals();
  
  d_diag_sef.reportName();
  d_diag_sef.reportFilename();
  d_diag_sef.setType( C_DBL );
  d_diag_sef.printType();
  d_diag_sef.setMode( DRA_RW );
  ddims.first=maxsef;
  ddims.second=1;
  d_diag_sef.setDimensions( ndims, ddims, rdims );
  if ( d_diag_sef.createDRA() != 0 )
    {
      cout << "failed to create the d_diag_sef: Aborting app " << endl;
      GA::error(" failed to create the d_diag_sef: Aborting app", GA::nodeid() );
    }
  d_diag_sef.printInternals();
}

/** Allocate the seven DRA data structures required to STORE GA to disk
    SUPPLEMENTAL version
 */
void AllocateDRA( HSOURCE state, int ndims, int maxsef, int maxsparse,int maxsparseBig, int maxwidth, PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, PsociDRArestart & d_cimat_supp,PsociDRArestart & d_icicol_supp,PsociDRArestart & d_number_supp )
{
  if ( state != CONSTRUCT ) {
    cout << "AllocateDRA-Supplemental form : apparent mismatch: state = " << state;
  }
  
  if (GA::nodeid() == 0 ) cout << "Preparing to Allocate DRAS for split GA-DRA Hamiltonian " << endl;
  
  pair<dra_size_t,dra_size_t> ddims;
  pair<dra_size_t,dra_size_t> rdims;
  
  rdims.first=-1;
  rdims.second=-1;
  
  d_cimat.reportName();
  d_cimat.reportFilename();
  d_cimat.setType( C_DBL );
  d_cimat.printType();
  d_cimat.setMode( DRA_RW );
  ddims.first=maxsef;
  ddims.second=maxsparse;
  d_cimat.setDimensions( ndims, ddims, rdims );
  if ( d_cimat.createDRA() != 0 )
    {
      cout << "failed to create the d_cimat: Aborting app " << endl;
      GA::error(" failed to create the d_cimat: Aborting app", GA::nodeid() );
    }
  d_cimat.printInternals();
  
  d_icicol.reportName();
  d_icicol.reportFilename();
  d_icicol.setType( C_INT );
  d_icicol.printType();
  d_icicol.setMode( DRA_RW );
  d_icicol.setDimensions( ndims, ddims, rdims );
  if ( d_icicol.createDRA() != 0 )
    {
      cout << "failed to create the d_icicol: Aborting app " << endl;
      GA::error(" failed to create the d_icicol: Aborting app", GA::nodeid() );
    }
  d_icicol.printInternals();
  
  d_number.reportName();
  d_number.reportFilename();
  d_number.setType( C_INT );
  d_number.printType();
  d_number.setMode( DRA_RW );
  ddims.first=maxsef;
  ddims.second=1;
  d_number.setDimensions( ndims, ddims, rdims );
  if ( d_number.createDRA() != 0 )
    {
      cout << "failed to create the d_number: Aborting app " << endl;
      GA::error(" failed to create the d_number: Aborting app", GA::nodeid() );
    }
  d_number.printInternals();
  
  d_diag_sef.reportName();
  d_diag_sef.reportFilename();
  d_diag_sef.setType( C_DBL );
  d_diag_sef.printType();
  d_diag_sef.setMode( DRA_RW );
  ddims.first=maxsef;
  ddims.second=1;
  d_diag_sef.setDimensions( ndims, ddims, rdims );
  if ( d_diag_sef.createDRA() != 0 )
    {
      cout << "failed to create the d_diag_sef: Aborting app " << endl;
      GA::error(" failed to create the d_diag_sef: Aborting app", GA::nodeid() );
    }
  d_diag_sef.printInternals();

/* Start opening up Supplemental (split) storage
*/

  maxsef = maxwidth;
  maxsparse = maxsparseBig;

  rdims.first=-1;
  rdims.second=-1;

  d_cimat_supp.reportName();
  d_cimat_supp.reportFilename();
  d_cimat_supp.setType( C_DBL );
  d_cimat_supp.printType();
  d_cimat_supp.setMode( DRA_RW );
  ddims.first=maxwidth;
  ddims.second=maxsparseBig;
  d_cimat_supp.setDimensions( ndims, ddims, rdims );
  if ( d_cimat_supp.createDRA() != 0 )
    {   
      cout << "failed to create the d_cimat_supp: Aborting app " << endl;
      GA::error(" failed to create the d_cimat_supp: Aborting app", GA::nodeid() );
    }   
  d_cimat_supp.printInternals();
  
  d_icicol_supp.reportName();
  d_icicol_supp.reportFilename();
  d_icicol_supp.setType( C_INT );
  d_icicol_supp.printType();
  d_icicol_supp.setMode( DRA_RW );
  d_icicol_supp.setDimensions( ndims, ddims, rdims );
  if ( d_icicol_supp.createDRA() != 0 ) 
    {   
      cout << "failed to create the d_icicol_supp: Aborting app " << endl;
      GA::error(" failed to create the d_icicol_supp: Aborting app", GA::nodeid() );
    }   
  d_icicol_supp.printInternals();
  
  d_number_supp.reportName();
  d_number_supp.reportFilename();
  d_number_supp.setType( C_INT );
  d_number_supp.printType();
  d_number_supp.setMode( DRA_RW );
  ddims.first=maxsef;
  ddims.second=1;
  d_number_supp.setDimensions( ndims, ddims, rdims );
  if ( d_number_supp.createDRA() != 0 ) 
    {   
      cout << "failed to create the d_number_supp: Aborting app " << endl;
      GA::error(" failed to create the d_number_supp: Aborting app", GA::nodeid() );
    }   
  d_number_supp.printInternals();
}

/** Specify and prepare to fetch the four DRA data structures from disk
 */
int SpecifyDRA( HSOURCE state, PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, int & readSparse)
{
  int retvalue=0;

  if ( state != RESTART ) {
    cout << "SpecifyDRA: apparent mismatch: state = " << state;
  }
  
  readSparse = 0;

  int ndims;
  pair<dra_size_t,dra_size_t> ddims;
  pair<dra_size_t,dra_size_t> rdims;

//CIMAT
  d_cimat.openDRAforREAD();
  d_cimat.inquireDRAsetREAD(); // ensures internal data structures are complete: needs to be refactored
  
  d_cimat.setMode( DRA_R );
  d_cimat.reportName();
  d_cimat.reportFilename();
  d_cimat.printType();
  d_cimat.printInternals();

  d_cimat.getDimensions( ndims, ddims, rdims );

  readSparse = ddims.second; // 

//ICICOL
  d_icicol.setMode( DRA_R );
  d_icicol.openDRAforREAD();
  d_icicol.inquireDRAsetREAD();
  d_icicol.getDimensions( ndims, ddims, rdims );
  
  d_icicol.reportName(); 
  d_icicol.reportFilename();
  d_icicol.printType();
  d_icicol.printInternals(); 
  
//NUMBER
  d_number.setMode( DRA_R );
  d_number.openDRAforREAD();
  d_number.inquireDRAsetREAD();
  d_number.getDimensions( ndims, ddims, rdims );

  d_number.reportName();
  d_number.reportFilename();
  d_number.printType();
  d_number.printInternals();

// DIAG_SEF
  d_diag_sef.setMode( DRA_R );
  d_diag_sef.openDRAforREAD();
  d_diag_sef.inquireDRAsetREAD();
  d_diag_sef.getDimensions( ndims, ddims, rdims );

  d_diag_sef.reportName();
  d_diag_sef.reportFilename();
  d_diag_sef.printType();
  d_diag_sef.printInternals();

  return( retvalue );
} 


/** Specify and prepare to fetch the four DRA data structures from disk
    Split version
 */
int SpecifyDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, 
                        PsociDRArestart & d_cimat_supp,PsociDRArestart & d_icicol_supp,PsociDRArestart & d_number_supp, int & readSparse, int & readSparse_supp, 
                        int & readWidth )
{
  int retvalue = 0;
  if ( state != RESTART ) {
    cout << "SpecifyDRA SPLIT: apparent mismatch: state = " << state;
  }
  
  readSparse = 0;
  readSparse_supp = 0;

  int ndims;
  pair<dra_size_t,dra_size_t> ddims;
  pair<dra_size_t,dra_size_t> rdims;

//CIMAT
  d_cimat.openDRAforREAD();
  d_cimat.inquireDRAsetREAD(); // ensures internal data structures are complete: needs to be refactored
  
  d_cimat.setMode( DRA_R );
  d_cimat.reportName();
  d_cimat.reportFilename();
  d_cimat.printType();
  d_cimat.printInternals();

  d_cimat.getDimensions( ndims, ddims, rdims );

  readSparse = ddims.second; // 


//ICICOL
  d_icicol.setMode( DRA_R );
  d_icicol.openDRAforREAD();
  d_icicol.inquireDRAsetREAD();
  d_icicol.getDimensions( ndims, ddims, rdims );
  
  d_icicol.reportName(); 
  d_icicol.reportFilename();
  d_icicol.printType();
  d_icicol.printInternals(); 
  
//NUMBER
  d_number.setMode( DRA_R );
  d_number.openDRAforREAD();
  d_number.inquireDRAsetREAD();
  d_number.getDimensions( ndims, ddims, rdims );

  d_number.reportName();
  d_number.reportFilename();
  d_number.printType();
  d_number.printInternals();

// DIAG_SEF
  d_diag_sef.setMode( DRA_R );
  d_diag_sef.openDRAforREAD();
  d_diag_sef.inquireDRAsetREAD();
  d_diag_sef.getDimensions( ndims, ddims, rdims );

  d_diag_sef.reportName();
  d_diag_sef.reportFilename();
  d_diag_sef.printType();
  d_diag_sef.printInternals();

//CIMAT_SUPP
  d_cimat_supp.openDRAforREAD();
  d_cimat_supp.inquireDRAsetREAD(); // ensures internal data structures are complete: needs to be refactored
  
  d_cimat_supp.setMode( DRA_R );
  d_cimat_supp.reportName();
  d_cimat_supp.reportFilename();
  d_cimat_supp.printType();
  d_cimat_supp.printInternals();

  d_cimat_supp.getDimensions( ndims, ddims, rdims );

  readWidth = ddims.first;
  readSparse_supp = ddims.second; //  

//ICICOL_SUPP
  d_icicol_supp.setMode( DRA_R );
  d_icicol_supp.openDRAforREAD();
  d_icicol_supp.inquireDRAsetREAD();
  d_icicol_supp.getDimensions( ndims, ddims, rdims );
  
  d_icicol_supp.reportName(); 
  d_icicol_supp.reportFilename();
  d_icicol_supp.printType();
  d_icicol_supp.printInternals(); 
  
//NUMBER_SUPP
  d_number_supp.setMode( DRA_R );
  d_number_supp.openDRAforREAD();
  d_number_supp.inquireDRAsetREAD();
  d_number_supp.getDimensions( ndims, ddims, rdims );

  d_number_supp.reportName();
  d_number_supp.reportFilename();
  d_number_supp.printType();
  d_number_supp.printInternals();

  return( retvalue );
} 

/** Send GA timer wrapped
 */
void SendHamiltonianToDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,
			   PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, 
			   PsociGAhamiltonian & hamilton, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  SendHamiltonianToDRA( state, d_cimat, d_icicol, d_number, d_diag_sef, hamilton );
  info.second = psociTime() - timein;
  info.first = g_rank;
}

void SendHamiltonianToDRA( HSOURCE state,  PsociDRArestart & d_cimat, PsociDRArestart & d_icicol,
                           PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, 
                           PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp,
                           PsociGAhamiltonian & hamilton, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  SendHamiltonianToDRA( state, d_cimat, d_icicol, d_number, d_diag_sef, d_cimat_supp, d_icicol_supp, d_number_supp, hamilton );
  info.second = psociTime() - timein;
  info.first = g_rank;
} 


/** Send GA data to disk
 */
void SendHamiltonianToDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,PsociDRArestart & d_number,PsociDRArestart & d_diag_sef, PsociGAhamiltonian & hamilton)
{
  if ( state != CONSTRUCT ) {
    cout << "mismatched call to SendHamiltonianToDRA? " << endl;
  }
  
  if ( GA::nodeid() == 0 ) cout << " Move GA to DRA " << endl;
  
  if ( d_cimat.dumpGAtoDRA( hamilton.fetchCIMAT() ) != 0 ) {
    cout << "Wholesale Dump CI mat failed " << endl;
  }
  d_cimat.flickDRA();
  if ( d_icicol.dumpGAtoDRA( hamilton.fetchICICOL() ) != 0 ) {
    cout << "Wholesale Dump ICICOL mat failed " << endl;
  }
  d_icicol.flickDRA();
  if ( d_number.dumpGAtoDRA( hamilton.fetchNUMBER() ) != 0 ) {
    cout << "Wholesale Dump NUMBER mat failed " << endl;
  }
  d_number.flickDRA();
  if ( d_diag_sef.dumpGAtoDRA( hamilton.fetchDIAG_SEF() ) != 0 ) {
    cout << "Wholesale Dump DIAG_SEF mat failed " << endl;
  }
  d_diag_sef.flickDRA();
  
  d_cimat.waitDRA();
  d_icicol.waitDRA();
  d_number.waitDRA();
  d_diag_sef.waitDRA();
  
}
/** Send GA data to disk SUPPLEMENTAL VERSION
 */
void SendHamiltonianToDRA( HSOURCE state, PsociDRArestart & d_cimat, PsociDRArestart & d_icicol, PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, PsociDRArestart & d_cimat_supp,
                          PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp, PsociGAhamiltonian & hamilton )
{
  if ( state != CONSTRUCT ) {
    cout << "mismatched call to SendHamiltonianToDRA-split? " << endl;
  }
  
  if ( GA::nodeid() == 0 ) cout << " Move Split GA to Split DRA " << endl;
  
  if ( d_cimat.dumpGAtoDRA( hamilton.fetchCIMAT() ) != 0 ) {
    cout << "Wholesale Dump CI mat failed " << endl;
  }
  d_cimat.flickDRA();
  if ( d_icicol.dumpGAtoDRA( hamilton.fetchICICOL() ) != 0 ) {
    cout << "Wholesale Dump ICICOL mat failed " << endl;
  }
  d_icicol.flickDRA();
  if ( d_number.dumpGAtoDRA( hamilton.fetchNUMBER() ) != 0 ) {
    cout << "Wholesale Dump NUMBER mat failed " << endl;
  }
  d_number.flickDRA();
  if ( d_diag_sef.dumpGAtoDRA( hamilton.fetchDIAG_SEF() ) != 0 ) {
    cout << "Wholesale Dump DIAG_SEF mat failed " << endl;
  }
  d_diag_sef.flickDRA();
  
  d_cimat.waitDRA();
  d_icicol.waitDRA();
  d_number.waitDRA();
  d_diag_sef.waitDRA();

  if ( d_cimat_supp.dumpGAtoDRA( hamilton.fetchCIMAT_supp() ) != 0 ) { 
    cout << "Wholesale Dump CI_supp mat failed " << endl;
  }
  d_cimat.flickDRA();
  if ( d_icicol_supp.dumpGAtoDRA( hamilton.fetchICICOL_supp() ) != 0 ) { 
    cout << "Wholesale Dump ICICOL_supp mat failed " << endl;
  }
  d_icicol.flickDRA();
  if ( d_number_supp.dumpGAtoDRA( hamilton.fetchNUMBER_supp() ) != 0 ) { 
    cout << "Wholesale Dump NUMBER_supp mat failed " << endl;
  }

  d_cimat_supp.waitDRA();
  d_icicol_supp.waitDRA();
  d_number_supp.waitDRA();
  
}

void FetchHamiltonianFromDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,
			      PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, 
			      PsociGAhamiltonian & hamilton, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  FetchHamiltonianFromDRA( state, d_cimat, d_icicol, d_number, d_diag_sef, hamilton );
  info.second = psociTime() - timein;
  info.first = g_rank;
}

void FetchHamiltonianFromDRA( HSOURCE state,  PsociDRArestart & d_cimat,PsociDRArestart & d_icicol,
                              PsociDRArestart & d_number, PsociDRArestart & d_diag_sef,
                              PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp,
                              PsociGAhamiltonian & hamilton, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  FetchHamiltonianFromDRA( state, d_cimat, d_icicol, d_number, d_diag_sef, d_cimat_supp, d_icicol_supp, d_number_supp, hamilton );
  info.second = psociTime() - timein;
  info.first = g_rank;
}


/** Fetch GA data to disk
 */

void FetchHamiltonianFromDRA( HSOURCE state, PsociDRArestart & d_cimat, PsociDRArestart & d_icicol, PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, PsociGAhamiltonian & hamilton)
{
  if ( state != RESTART ) {
    cerr << "Erroneous call to FetchHamiltonianToDRA " << endl;
    GA::error( "Erroneous call to FetchHamiltonianToDRA ", state);
  }
  
  if (GA::nodeid() == 0 ) cout << " Move DRA to GA " << endl;
  
  if ( d_cimat.readGAfromDRA( hamilton.fetchCIMAT())  != 0 ) {
    cerr << "Failed to read CIMAT " << endl;
    GA::error("Failed to read CIMAT ", -1 );
  }
  if ( d_icicol.readGAfromDRA( hamilton.fetchICICOL())  != 0 ) {
    cerr << "Failed to read ICICOL " << endl;
    GA::error("Failed to read ICICOL ", -1 );
  }
  if ( d_number.readGAfromDRA( hamilton.fetchNUMBER())  != 0 ) {
    cerr << "Failed to read NUMBER " << endl;
    GA::error("Failed to read NUMBER ", -1 );
  }
  if ( d_diag_sef.readGAfromDRA( hamilton.fetchDIAG_SEF())  != 0 ) {
    cerr << "Failed to read DIAG_SEF " << endl;
    GA::error("Failed to read DIAG_SEF ", -1 );
  }
  
  d_cimat.waitDRA();
  d_icicol.waitDRA();
  d_number.waitDRA();
  d_diag_sef.waitDRA();
}

/* Split DRA-GA fetch method
*/
void FetchHamiltonianFromDRA( HSOURCE state, PsociDRArestart & d_cimat, PsociDRArestart & d_icicol, PsociDRArestart & d_number, PsociDRArestart & d_diag_sef, 
        PsociDRArestart & d_cimat_supp, PsociDRArestart & d_icicol_supp, PsociDRArestart & d_number_supp, PsociGAhamiltonian & hamilton)
{
  if ( state != RESTART ) {
    cerr << "Erroneous call to Split version FetchHamiltonianToDRA " << endl;
    GA::error( "Erroneous call to Split version of FetchHamiltonianToDRA ", state);
  }
  
  if (GA::nodeid() == 0 ) cout << " Move DRA to GA " << endl;
  
  if ( d_cimat.readGAfromDRA( hamilton.fetchCIMAT())  != 0 ) {
    cerr << "Failed to read CIMAT " << endl;
    GA::error("Failed to read CIMAT ", -1 );
  }
  if ( d_icicol.readGAfromDRA( hamilton.fetchICICOL())  != 0 ) {
    cerr << "Failed to read ICICOL " << endl;
    GA::error("Failed to read ICICOL ", -1 );
  }
  if ( d_number.readGAfromDRA( hamilton.fetchNUMBER())  != 0 ) {
    cerr << "Failed to read NUMBER " << endl;
    GA::error("Failed to read NUMBER ", -1 );
  }
  if ( d_diag_sef.readGAfromDRA( hamilton.fetchDIAG_SEF())  != 0 ) {
    cerr << "Failed to read DIAG_SEF " << endl;
    GA::error("Failed to read DIAG_SEF ", -1 );
  }
  
  d_cimat.waitDRA();
  d_icicol.waitDRA();
  d_number.waitDRA();
  d_diag_sef.waitDRA();

  if ( d_cimat_supp.readGAfromDRA( hamilton.fetchCIMAT_supp())  != 0 ) {
    cerr << "Failed to read CIMAT_supp " << endl;
    GA::error("Failed to read CIMAT_supp ", -1 );
  }
  if ( d_icicol_supp.readGAfromDRA( hamilton.fetchICICOL_supp())  != 0 ) {
    cerr << "Failed to read ICICOL_supp " << endl;
    GA::error("Failed to read ICICOL_supp ", -1 );
  }
  if ( d_number_supp.readGAfromDRA( hamilton.fetchNUMBER_supp())  != 0 ) {
    cerr << "Failed to read NUMBER_supp " << endl;
    GA::error("Failed to read NUMBER_supp ", -1 );
  }

  d_cimat_supp.waitDRA();
  d_icicol_supp.waitDRA();
  d_number_supp.waitDRA();


}








