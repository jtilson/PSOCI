/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: A rudementary class for managing servicess for DRA: Basically interfaces 
              to the DRA library 

 History:

 

**************************************************************************************/
/**
 *   @file PsociDRAservices.C
 *
 */

#include "PsociTimer.hpp"
#include "PsociDRAservices.hpp"

// Class definitions

using namespace std;

PsociDRAservices::PsociDRAservices()
{
  PSOCI_MAX_ARRAYS = DEF_PSOCI_MAX_ARRAYS; 
  PSOCI_MAX_ARRAYSIZE = DEF_PSOCI_MAX_ARRAYSIZE;
  PSOCI_TOTALDISKSPACE = DEF_PSOCI_TOTALDISKSPACE;
  PSOCI_MAX_MEMORY = DEF_PSOCI_MAX_MEMORY;
}

// Collective
void PsociDRAservices::exitDRA()
{
  if ( GA::nodeid() == 0 ) cout << "Terminating DRA services" << endl;
  if ( DRA_Terminate() != 0 ) {
    cerr << "Unsuccessful termination of the DRA layer: Aborting job " << endl;
    GA::Terminate();
  }
}

//Collective call  Init DRA values
void PsociDRAservices::initDRA( int rank )
{
  DRA_status = DRA_Init( PSOCI_MAX_ARRAYS, PSOCI_MAX_ARRAYSIZE, PSOCI_TOTALDISKSPACE,
            PSOCI_MAX_MEMORY );
  //cout << "DRA services status is " << DRA_status << endl;
  PSOCI_NUM_DRA_IO_FILES = GA::nodes();
  PSOCI_NUM_DRA_IO_SERVERS = GA::nodes();
}

// DRA doesn't permit numfiles <> numservers
//Collective call  Init DRA values
void PsociDRAservices::initDRA( int rank, int servers )
{
  DRA_status = DRA_Init( PSOCI_MAX_ARRAYS, PSOCI_MAX_ARRAYSIZE, PSOCI_TOTALDISKSPACE,
            PSOCI_MAX_MEMORY );
  //cout << "DRA services status is " << DRA_status << endl;
  // PSOCI_NUM_DRA_IO_FILES = GA::nodes();
  PSOCI_NUM_DRA_IO_FILES = servers;
  PSOCI_NUM_DRA_IO_SERVERS = servers;
  //cout << "vals are " << PSOCI_NUM_DRA_IO_FILES<<" "<<PSOCI_NUM_DRA_IO_SERVERS<<endl;
  DRA_Set_default_config( PSOCI_NUM_DRA_IO_FILES, PSOCI_NUM_DRA_IO_SERVERS );
// How are thes eprocesses distributed about nodes?
}


//Local call
void PsociDRAservices::reportDRAIOservers( void )
{
   if ( GA::nodeid() == 0 ) {
      cout << " PSOCI_NUM_DRA_IO_FILES = " << PSOCI_NUM_DRA_IO_FILES;
      cout << " PSOCI_NUM_DRA_IO_SERVERS = " << PSOCI_NUM_DRA_IO_SERVERS << endl << endl;
   }
}


//Local call
void PsociDRAservices::reportDRAParams()
{
    if ( GA::nodeid() == 0 ) {
    cout << " Current settings for DRA_init " 
         << " PSOCI_MAX_ARRAYS = " << PSOCI_MAX_ARRAYS
         << " PSOCI_MAX_ARRAYSIZE(B) = " << PSOCI_MAX_ARRAYSIZE
         << " PSOCI_TOTALDISKSPACE(B) = " << PSOCI_TOTALDISKSPACE
         << " PSOCI_MAX_MEMORY(B) = " << PSOCI_MAX_MEMORY << endl << endl;
    }
}

