/**************************************************************************************

* Copyright (c) 2010,2012 RENCI.
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
 *   @file PsociTaskManager.hpp
 *
 */

#ifndef PSOCITASKMANAGER_H
#define PSOCITASKMANAGER_H

#ifdef TCGMSG 
#include <tcgmsg.h>
#endif

#include "PsociTimer.hpp"
#include "PsociDeterminants.hpp"
#include "PsociHamiltonian.hpp"
#include "PsociGAhamiltonian.hpp"
#include "PsociBlasLapack.hpp"
#include "PsociTaskManager.hpp"

using namespace std;

class PsociTaskManager {
  
private:

  int l_nxtval_chunk;
  GA::GlobalArray * g_nxtval;

  int nleft;
  int icount;

public:

  explicit PsociTaskManager(int chunk, GA::GlobalArray * main_g_nxtval);
  int fetchNxtvalChunk(void );
  void setNxtvalChunk( int nchunk );
  long nxtask( int nproc ); 
  long nxtaskRev( int nproc, int topvalue );
  long nxtaskAlt( int nproc );
  void initNxtask();
  void initNxtask( int topvalue );
  void destroyNxtask();
  void nxtaskReset();
  void nxtaskReset( int topval );


};

#endif //PSOCITASKMANAGER_H
