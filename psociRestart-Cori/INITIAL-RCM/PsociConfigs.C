
//TODO remove the superfluous GA::sync() calls used for debugging
/**************************************************************************************

* Copyright (c) 2010,2011 RENCI.
* All rights reserved. This program and the accompanying materials
* MAY BE available under the terms of the RENCI Open Source License
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license


 Classes: 

 Description: Read process and validate spatial configs used for the SOCI
              For now we assume PARALLEL access to these even though most of the data
              are handled by one node.

 History:


**************************************************************************************/
/**
 *   @file PsociConfigs.C
 *
 */


/* TODO scramble the order of the configurations for better load balance
   When nxtval_chunk is big, we get poor distribution; perhaps because of the
   decidedly  no-random  configuration order ( all jsym=1,all jsym=2,etc.)

   This probably requires a new pack generate routimner that uses a Vector instead of a multimap
   As long as the "index" values are unique and correspond to the GA index, order doesn't matter
*/


#include "PsociTimer.hpp"
#include "PsociConfigs.hpp"

// Class definitions

using namespace std;

PsociConfigs::PsociConfigs()
{
  local_filename = DEFAULT_FILENAME;
}

PsociConfigs::PsociConfigs( string filename )
{
  local_filename = filename;
}

int  PsociConfigs::fetchSpinParity()
{
    return( spin_parity );
}

int  PsociConfigs::fetchMaxOpenShells()
{
    return( n_openshells );
}

int  PsociConfigs::fetchStateSymmetry()
{
    return( st_symmetry );
}

int  PsociConfigs::fetchNumBasisFtns()
{
    return( n_basisftns );
}

int  PsociConfigs::fetchNumElectrons()
{
    return( n_electrons );
}


//Destroy the current multimap: not needed once the GA array is built
void PsociConfigs::tearDownConfigs()
{
  if (GA::nodeid() == 0 ) cout << "tearing down configs " << endl;
  configs.clear();
  vec_configs.clear();
  vec_configs_index.clear();
  n_configurations = 0;
#ifdef SINGLENODECONFIGS 
  destroyGAconfigs();
#endif
}

void PsociConfigs::printFilename()
{
  cout << "Current configuration filename is " << local_filename << endl;
}

void PsociConfigs::OpenFile()
{
  ilocalfile.open( local_filename.c_str(), ios::in ); 
  if (!ilocalfile.is_open()) {
    GA::error("Cannot open PSOCI configuration input file ",1);    
    cerr << "Cannot open PSOCI configuration input file " << local_filename << " Aborting read " << endl;
    exit(1);
  }
}

//
void PsociConfigs::CloseFile()
{
  ilocalfile.close();
}

//Timer wrapped version
int PsociConfigs::readConfigs( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int status = readConfigs();
  info.second =  psociTime() - timein;
  info.first = g_rank;
  return( status );
}
    
//Return 0 on success else returns input line number that error occured or simply -1;
int PsociConfigs::readConfigs()
{
  configs.clear();
  set< pair<int,string> > temp_set; // used to remove redundancies. NOTE same configs with diff syms are not caught.
  pair< set<pair<int,string> >::iterator, bool> pcheck;
  bool rem_dup=false;

  if ( GA::nodeid() != 0 ) return(0);
  OpenFile();
  if (!ilocalfile.is_open()) {
    GA::error("Cannot open CONFIGURATION input file ",1);
    return(-1);
  }
  const int MAXLINE=256;
  char ttl[MAXLINE];
  
  // Get job title
  
  ilocalfile.getline(ttl, MAXLINE);
  if ( !ilocalfile ) {
    cout << "Failed to get a proper title " << ttl << endl;
    CloseFile();
    return(1);
  } else {
    string title( ttl );
    local_title = title;
  }
  
  // Read the job parmaters
  
  char line[MAXLINE];
  int nbf, nelec, nsym, ncfgs, nopen;
  ilocalfile.getline(line, MAXLINE);
  stringstream wss(line);
  if ( wss >> nbf >> nelec >> nsym >> ncfgs >> nopen ) {
    n_openshells = nopen;
    n_basisftns = nbf;
    n_electrons= nelec;
    st_symmetry = nsym;
    n_configurations = ncfgs;
    spin_parity = ( n_electrons % 2 ); //Compute now. No need to check this again
  } else {
    GA::error("Failed to read 5 ints specifying the job: return as a failure ",1);
    CloseFile();
    return(2);
  }
  
  // Read the configs - perform some basic validation as well ( but not symmetry that uses its own method )
  int isym;
  int iopen;
  int ielec;
  string iconf;
  int numcfgs=0;
  
  int temp_openshells = 0;
  
  while ( ilocalfile.getline(line, MAXLINE) ) {
    ++numcfgs;
    stringstream cfg(line);
    if ( !(cfg >> isym) ) {
      GA::error("Failed to read a symmetry specification for a configuration: ",1);
      CloseFile();
      return(numcfgs+2);
    }
    if ( !(cfg >> iconf) ) {
      GA::error("Failed to read a configuration at numcfgs = ",1);
      CloseFile();
      return(numcfgs+2);
    }
    
    // Check size of the config
    
    if ( iconf.size() != n_basisftns ) {
      GA::error("Error: Found wrong config size. Expected n_basisftns = ",1);
      cerr << " Found a configs of size " << iconf.size() << " at numcfgs = " <<  numcfgs;
      CloseFile();
      return( numcfgs+2 ); 
    }
    
    // check current number of openshells 
    iopen=0;
    for (int i=0; i<iconf.size(); ++i) {
      if (iconf[i] == '1' ) ++iopen;
    }

    if (iopen > n_openshells ) {
      GA::error("Error: Found more open shells than expected: Expected  maximum of ",n_openshells);
      cerr << " at config number " << numcfgs << endl;
      CloseFile();
      return( numcfgs+2 );
    }
    temp_openshells = max( temp_openshells, iopen );

    // Might as well check the number of elecs as well
    
    ielec=0;
    for (int i=0; i<iconf.size(); ++i) {
      if (iconf[i] == '1' ) ++ielec;
      if (iconf[i] == '2' ) ielec += 2;
      if (iconf[i] != '0' && iconf[i] != '1' && iconf[i] != '2' ) {
         GA::error("Found on orbital specification not 0,1, or 2: at cfg number ",numcfgs);
         cerr << " Aborting: " << endl;
         GA::error("Found on orbital specification not 0,1, or 2: at cfg number ", i);
      }
    }

    if ( iopen == 11 ) {
       for(int yy=0; yy<n_basisftns ; ++yy) {
       cout << iconf[yy];
    }
    cout << isym<<endl;
    }

    if (ielec != n_electrons ) {
      GA::error("Error: Found different electrons than expected: found ",ielec);
      cerr << " at config number " << numcfgs << endl;
      CloseFile();
      return( numcfgs+2 );
    }
    
    if ( !(ielec-iopen) % 2 ) {
      GA::error(" ielec-iopen % 2 is not zero: ",1);
      return( numcfgs+2 );
    }

    //Okay store the current confuguration
    
//    configs.insert(pair<const int,string>(isym, iconf));
      pcheck = temp_set.insert( pair<const int,string>(isym, iconf)); //quietly removes dups
      if ( pcheck.second != true ) {
          cout << "Duplicate config not inserted " << isym <<" "<<iconf << endl;
          --numcfgs;
          rem_dup=true;
      }
  } 
  
  /* I didn't want to adjust UP as this may have other consequences: Adjusting down however
     should be fairly safe
  */

   set< pair<int,string> >::iterator sit;
   for(sit = temp_set.begin(); sit != temp_set.end(); ++sit) {
      configs.insert( *(sit) );
   }

  // Check maximum openshells again

  if ( temp_openshells < n_openshells ) {
    cout << "Processing determined the actual maximum num open shells is " <<  temp_openshells;
    cout << " Though the input file had provided a value of " << n_openshells;
    cout << " Adjusting DOWN the value of n_openshells " << endl;
    n_openshells = temp_openshells;
  }
  
  // Check number of configs
  
  if ( configs.size() != n_configurations ) {
    cout << "Incorrect number of configs read. Expected " << n_configurations << "Found ";
    cout << " numcfs = " <<  numcfgs << " Yielding a map of size " << configs.size()  << endl;
    n_configurations = configs.size();
    CloseFile();
/* In this case if dups were removed then its impossible to have predicted what n_configurations
   would be. SO accept the outcome. On the other hand if NO dups were found then the user
   probably simply mistyped the value and thus we should stop!
*/
    if (rem_dup == true) {
      cout << " duplicates have been found and removed: permitting continuation " << endl << endl;
      return( 0 );
    } else {
    return( numcfgs + 2 );
    }
  }
  CloseFile();
  return( 0 );
}

int PsociConfigs::numTotalSpatials()
{
  return( n_configurations );
}

int PsociConfigs::readConfigsToGA( pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  int status = readConfigsToGA();
  info.second =  psociTime() - timein;
  info.first = g_rank;
  return( status );
}

//Return 0 on success else returns input line number that error occured or simply -1;

/* A new method. This one starts with rank=0 doing the read and then immediately pushing
   contents to GA. Later we wil fiogure out how to parallel process them

   will still need to do a subsequent brdcstParams so everyone knows the basic parameters

   GA space MUST BE preallocate here as only rank=0 is allowed in .

   This is a horrible combination of parallel and rank=0 work
*/

//Return 0 on success else returns input line number that error occured or simply -1;
// THIS IS A PARALLEL method all nodes must respond

/* Currently this method CANNOT check for duplicates. User beware ! */

int PsociConfigs::readConfigsToGA()
{
  const int MAXLINE=256;
  char ttl[MAXLINE];
  int nbf, nelec, nsym, ncfgs, nopen;
  char line[MAXLINE];

  if ( GA::nodeid() == 0 ) {
  OpenFile();
  if (!ilocalfile.is_open()) {
    GA::error("readConfigsToGA: Cannot open CONFIGURATION input file",1 );
    return(-1);
  }
  
  // Get job title
  
  ilocalfile.getline(ttl, MAXLINE);
  if ( !ilocalfile ) {
    cout << "Failed to get a proper title " << ttl << endl;
    CloseFile();
    return(1);
  } else {
    string title( ttl );
    local_title = title;
  }
  
  // Read the job parmaters
  
  ilocalfile.getline(line, MAXLINE);
  stringstream wss(line);
  if ( wss >> nbf >> nelec >> nsym >> ncfgs >> nopen ) {
    n_openshells = nopen;
    n_basisftns = nbf;
    n_electrons= nelec;
    st_symmetry = nsym;
    n_configurations = ncfgs;
    spin_parity = ( n_electrons % 2 ); //Compute now. No need to check this again
  } else {
    cout << "Failed to read 5 ints specifying the job: return as a failure " << endl;
    CloseFile();
    return(2);
  }
  
  } // not rank=0
  
  /* At this point we can open up the GA space. then all codes <> 0 can leave 
     array is [86][numconfs] is size
  */
  
  // can't create the GA 'till we know how big.
  
  GA::brdcst( &nbf, sizeof(int), 0) ;
  GA::brdcst( &n_configurations, sizeof(int), 0) ;
  
  // cout << "create configs " << endl;
  createGAconfigs( nbf, n_configurations );
  
  if ( GA::nodeid() != 0 ) return(0); // not needed anymore. This needs serious refactoring.
  
  // Read the configs - perform some basic validation as well ( but not symmetry that uses its own method )
  int isym;
  int iopen;
  int ielec;
  string iconf;
  int numcfgs=0;
  
  int temp_openshells = 0;
  int indexGA=0;
  
  while ( ilocalfile.getline(line, MAXLINE) ) {
    ++numcfgs;
    stringstream cfg(line);
    if ( !(cfg >> isym) ) {
      cerr << "Failed to read a symmetry specification for a configuration: numcfgs = " << numcfgs << endl;
      CloseFile();
      return(numcfgs+2);
    }
    if ( !(cfg >> iconf) ) {
      cerr << "Failed to read a configuration at numcfgs = " << numcfgs << endl;
      CloseFile();
      return(numcfgs+2);
    }
    
    // Check size of the config
    
    if ( iconf.size() != n_basisftns ) {
      cerr << "Error: Found wrong config size. Expected n_basisftns = " << n_basisftns;
      cerr << " Found a configs of size " << iconf.size() << " at numcfgs = " <<  numcfgs;
      CloseFile();
      return( numcfgs+2 ); 
    }
    
    // check current number of openshells 
    iopen=0;
    for (int i=0; i<iconf.size(); ++i) {
      if (iconf[i] == '1' ) ++iopen;
    }
    if (iopen > n_openshells ) {
      cerr << "Error: Found more open shells than expected: Expected  maximum of " << n_openshells << " found " << iopen;
      cerr << " at config number " << numcfgs << endl;
      CloseFile();
      return( numcfgs+2 );
    }
    temp_openshells = max( temp_openshells, iopen );
    
    // Might as well check the number of elecs as well
    
    ielec=0;
    for (int i=0; i<iconf.size(); ++i) {
      if (iconf[i] == '1' ) ++ielec;
      if (iconf[i] == '2' ) ielec += 2;
      if (iconf[i] != '0' && iconf[i] != '1' && iconf[i] != '2' ) {
	cerr << "Found on orbital specification not 0,1, or 2: at cfg number " << numcfgs << " index " << i << endl;
	cerr << " Aborting: " << endl;
	GA::error("Found on orbital specification not 0,1, or 2: at cfg number ", i);
      }
    }
    
    if (ielec != n_electrons ) {
      cerr << "Error: Found different electrons than expected: Expected " << n_electrons << " found " << ielec;
      cerr << " at config number " << numcfgs << endl;
      CloseFile();
      return( numcfgs+2 );
    }
    
    if ( !(ielec-iopen) % 2 ) {
      cerr << " ielec-iopen % 2 is not zero: " << ielec << " " << iopen << endl;
      return( numcfgs+2 );
    }
    
    //Okay store the current confuguration
    
/*
    int strSize = iconf.size();
    vector<char> conticonf(0); // resolve intel11/12 problems with padding on a string
    for(int k=0; k<strSize ; ++k ) {
       conticonf.push_back( iconf[k] );
    }
*/

    vector<int> confgs;
    convertStringToIntArray(iconf,confgs);
    pushConfigsToGA( indexGA, nbf, isym, confgs ); 
    ++indexGA;
  }
  
  // Check maximum openshells again

  if ( temp_openshells < n_openshells ) {
    cout << "Processing determined the actual maximum num open shells is " <<  temp_openshells;
    cout << " Though the input file had provided a value of " << n_openshells;
    cout << " Adjusting DOWN the value of n_openshells " << endl;
    n_openshells = temp_openshells;
  }
  
  // Check number of configs
  
  CloseFile();
  return( 0 );
}

int PsociConfigs::numBasisFunctions()
{
  return( n_basisftns );
}

void PsociConfigs::printParams()
{
  cout << "    Number of basis functions = " << n_basisftns << endl;
  cout << "          Number of electrons = " << n_electrons << endl;
  cout << "Maximum number of open shells = " << n_openshells << endl;
  cout << "     Number of configurations = " << n_configurations << endl;
  cout << "              State symmetry is " << st_symmetry << endl;
  cout << "                 Spin Parity is " << spin_parity << endl;
  cout << endl;
}

/* in the order nbf, nelec, st_symmetry, ncfgs, nopen to
   conform with the order observed in the original input file */

void PsociConfigs::setParams( vector<int> & params )
{
  if ( params.size() < 6 ) {
    cout << "Warning: expected a params vector with 6 or more entries: returning " << params.size() << endl;
    GA::Terminate();
  }
  if ( params.size() != 6 ) {
    cout << "Warning: expected a params vector of size 6: Found one of size " << params.size() << endl;
  }
  n_basisftns      = params[0];
  n_electrons      = params[1];
  st_symmetry      = params[2];
  n_configurations = params[3];
  n_openshells     = params[4];
  spin_parity      = params[5];
}

//
void PsociConfigs::fetchParams( vector<int> & params )
{
  if ( params.size() != 0 ) {
    cout << "Warning: expected an empty params vector: clearing it out" << endl;
  }
  params.clear();
  params.push_back( n_basisftns );
  params.push_back( n_electrons );
  params.push_back( st_symmetry );
  params.push_back( n_configurations );
  params.push_back( n_openshells );
  params.push_back( spin_parity );
}

/* In the method fetchConfigAlt, we want to be able to access the "configs"
   object using the index because the loop strucvutre used ( like in fetchConfigs )
   will not scale well for the pseudo-direct codes. So since that array will not be too big
   we will simply use the API to copy multimap into vector.
   Use the available extra call to destroy the multimap if you need space 
*/

// NOTE we may need to prune memory usage here
void PsociConfigs::copyInternalMapToVector()
{
  int startsize = configs.size(); // configs could be zero if using SCRAMBLE methods
  int stopsize = 0;

#ifdef CONFIGSCRAMBLE
  // If using scramble then configs is now empty BUT VectorToVectormap already populated vec_configs
  // so just skip the work
  if ( GA::nodeid() == 0 ) cout << "copyInternalMapToVector:SCRAMBLE approach " << vec_configs.size()<< endl;
  
#else
  startsize = configs.size();
  vec_configs.clear();
  multimap<int,string>::iterator it;
  it = configs.begin();
  while ( it != configs.end() ) {
    vec_configs.push_back( (*it) );
    ++it;
  }
#endif

  stopsize = vec_configs.size();
  
  if ( GA::nodeid() == 0 ) {
    if ( stopsize != 0 ) {
      cout << "Multimap to Vector copy successful "<< stopsize<<endl;
    } else {
      GA::error("Multimap to Vector copy failed: map-vec is", startsize-stopsize);
    }
  }
}

/*Mostly for debugging
  This version is for use by replicated configs-based methods.
  Generally you only need one core to call but no constraints 
  are imposed
*/
void PsociConfigs::printConfigurations()
{
  cout << endl << "Requested to printout entire configuration list" << endl;
  cout << "Note: Configurations have been re-ordered by symmetry " << endl;
  if ( n_configurations <= MAX_CONFIG_PRINT ) { 
    cout << "     printConfigurations:" << endl;
    cout << "    -- isym  configuration -- " << endl;
    multimap<int,string>::iterator it, out, inner;
    
    it = configs.begin();
    while ( it != configs.end() )
      {
        int isym = (*it).first;
        out = configs.upper_bound(isym);
        for( inner=it; inner != out; ++inner )
          {
            cout << isym << " " << (*inner).second << " " << endl;
          }
        cout << endl;
        it = out;
      } 
  } else {
    cout << "Configs list too large to print: is of size " << n_configurations;
    cout << " But must be less than MAX_CONFIG_PRINT = " << MAX_CONFIG_PRINT  << endl;
  }
  return;
}

/* A print method for the new distributed configuration methods
   ALL cores must participate
*/

void PsociConfigs::printDistributedConfigurations()
{
#ifndef SINGLENODECONFIGS
  GA::error("Application must be compiled with -DSINGLENODECONFIGS to access this method: printDistributedConfigurations",-1);
#endif
  
  int g_rank = GA::nodeid();
  int g_size = GA::nodes();
  
  if ( GA::nodeid() == 0 ) {
    cout << endl << "Requested to printout entire (distributed) configuration list" << endl;
    cout << endl;
  }
  GA::sync();
  
  if ( n_configurations > MAX_CONFIG_PRINT && g_rank == 0 ) { 
    cout << "Configs list too large to print: is of size " << n_configurations;
    cout << " But must be less than MAX_CONFIG_PRINT = " << MAX_CONFIG_PRINT  << endl;
    return;
  }
  
  //Print configs out in order one core at a time
  
  vector<pair<int,string> >::iterator it;
  
  for( int i=0; i< g_size; ++i ) {
    if ( i == g_rank ) {  
      if ( i == 0 ) {
	cout << "     printConfigurations" << endl;
	cout << "    -- isym  configuration -- " << endl;
      }
      
      if ( n_configurations <= MAX_CONFIG_PRINT ) {
	cout << "local configs: my rank is " << g_rank << endl;
	
	for( it=vec_configs.begin(); it!= vec_configs.end(); ++it) 
	  {
	    cout << (*it).first <<" "<< (*it).second << endl;
	  } 
	GA::sync();
      }
    }
  } 
  return;
}

//Timer wrapped brdcst
void PsociConfigs::brdcstConfigs( int root, pair<int,double> & info )
{
  int g_rank = GA::nodeid();
  double timein = psociTime();
  brdcstConfigs( root );
  info.second =  psociTime() - timein;
  info.first = g_rank;
}

/* We need a better data structure for handling the configs in a parallel environment
   Namely, the need to easuly Brdcast the data to/fro. Multimaps have no guarentee of
   being contiguous. So we should simply convert to vectors which ARE. Note, 'pair' is
   not contiguous either so we cannot use that.
*/

// For now create two arrays: one for isyms and one for configs. We also want to know how many bytes we are using ?
// Do we simply assume all configuration strings are the same size?

// We anticipate rank=0 doing the work, converting to vectors brdcasting to all and all 
// back transforming them to multimaps.

// Serialization step: Need something more specific controlled for the BRDCST step

//Set this up as a COLLECTIVE operation with the serialization/deserialization internal. 
//This should be called be all nodes AFTER a read configs call

//root_node must refer to the ONE node who actually READ the configs file

void PsociConfigs::brdcstConfigs( int root )
{
  int g_rank = GA::nodeid();
  if ( root > GA::nodes() ) {
    cerr << "brdcstConfig: Failure: root > g_size: root is " << root << " size is " << GA::nodes() << endl;
    GA::Terminate();
  }
  
// Grab currently set params on root and brdcst out to all

  vector<int> params(6);
  if( g_rank == root ) { 
    fetchParams( params );
  }
  GA::brdcst( &params[0], sizeof(int)*params.size(), root );
  setParams( params );
  
  vector<int> isym;
  vector<short> confs;
  vector<long> sizes(2); // refers to confs and isym arrays
  
  // Have root generate the serialized objects
  // and distribute the sizes of the serialized data objects (remember we converted to shorts)

  if ( g_rank == root ) {
    convertMultimapToVector( isym, confs, sizes );
  }
  GA::brdcst( &sizes[0], sizeof(long)*2, 0) ;
  
#if 0
   cout << g_rank << "sym in elements is " << sizes[0]/sizeof(int) <<  " sym in Bytes " << sizes[0] << endl;
   cout << g_rank << "confs in elements is " << sizes[1]/(nbasis*sizeof(short)) << " confs in Bytes " << sizes[1] << endl;
   cout << g_rank << " Attempt to BRDCST the serialized data to all nodes " << endl;
#endif

   // Have all REMOTE guys allocate space and BRDCST out the serialized objects

   if ( g_rank != root ) { //Root is already allocated
     isym.resize( sizes[0]/sizeof(int) );
     confs.resize( sizes[1]/sizeof(short) );
   }

#ifdef CONFIGSCRAMBLE
// Shuffle the contents of isym and confs. converttoVectormap will retain order

  if ( GA::nodeid() == root ) {
  shuffleConfigurations(root, isym, confs );
  }

#endif

   GA::brdcst( &isym[0], sizes[0], root );
   GA::brdcst( &confs[0], sizes[1], root );
   
#ifdef CONFIGSCRAMBLE
   // At this stage everyone has the replicated configutation lists.

     convertVectorToVectormap( isym, confs ); // populates vec_configs
   if ( g_rank == root ) cout << "Using scramble code: deleting multimap object" << endl;
   configs.clear(); // Empty multimap data 

#else
   // Lastly everyone needs to deserialize their data back into the maps
   if ( g_rank != root ) { // local multimap already populated
     convertVectorToMultimap( isym, confs );
   }
#endif
}

void PsociConfigs::brdcstParams( int root )
{
  int g_rank = GA::nodeid();
  if ( root > GA::nodes() ) {
    cerr << "brdcstParams: Failure: root > g_size: root is " << root << " size is " << GA::nodes() << endl;
    GA::Terminate();
  }
  
  // Grab currently set params on root and brdcst out to all
  
  vector<int> params(6);
  if( g_rank == root ) { 
    fetchParams( params );
  }
  GA::brdcst( &params[0], sizeof(int)*params.size(), root );
  setParams( params );
  
}

/** Uhg lots of problems handling the strings for MPI/GA brdcst. Since this only happens once and since
    the data structure is temporary, for serialization convert string to an array of SHORTS and back again

    Further all orbitals have already been checked to be (0, 1, or 2) at the read and so this is safe
    each vector<short> will be sizeof(short)*nbf long
*/

void PsociConfigs::convertMultimapToVector( vector<int> & isym, vector<short> & confs, vector<long> & sizes )
{
  if ( isym.size () >= 1 || confs.size() >= 1) {
    cout << "Warning: convertMaptoVector: expected empty vector objects. Will empty them now " << endl;
  }
  isym.clear();
  confs.clear();
  sizes.clear();
  
  long isymSizeB=0;
  long confSizeB=0;
  
  multimap<int,string>::iterator it, out, inner;
  
  it = configs.begin();
  while ( it != configs.end() )
    {
      int sym = (*it).first;
      out = configs.upper_bound(sym);
      for( inner=it; inner != out; ++inner )
	{
	  isymSizeB += sizeof( sym );
	  isym.push_back( sym );
	  
          //vector<short> vshort;
          convertStringToShortArray( (*inner).second, confs );
	  
	  //confSizeB += vshort.size()*sizeof(short); 
	  //confs.push_back( vshort );  not neede any more
	}
      cout << endl;
      it = out;
    }
  confSizeB += confs.size()*sizeof(short);
  sizes.push_back( isymSizeB );
  sizes.push_back( confSizeB );
#if 0
  cout << "Conversion to Vector: Size (B) of sym is " << isymSizeB << " size (B) of linearized conf is " << confSizeB << endl;
  cout << "Conversion to Vector: Word in sym is " << isym.size() << "   words in confs is " << confs.size()/nbasis << endl;
#endif
}

//Conversion of presumably MPI sent vector of data back into a local multimap array

// Deserialization step
void PsociConfigs::convertVectorToMultimap( vector<int> & isym, vector<short> & confs )
{
  int nbasis = n_basisftns;
  if ( nbasis == 0 ) { 
    cerr << "convertVectorToMultimap: nbasis is zero cannot recover " << endl;
    GA::Terminate();
  }
  if ( isym.size() == 0 || confs.size() == 0) {
    cout << "Warning: convertVectorToMap: expected non-empty vector objects. Need to abort" << endl;
    GA::Terminate();
  }
  //Clear out the existing configs multimap;
  configs.clear();
  
  
  // Are confs and isym the same length?
  if ( confs.size()/nbasis != isym.size() ) {
    cerr << nbasis << " convertVectorToMap: Non comforming confs and sym arrays: confs (estimated) is of size " << confs.size()/nbasis;
    cerr << " and sym is of size " << isym.size() << " Cannot recover " << endl;
    return;
  }
  if ( confs.size()/nbasis != isym.size() ) {
    cerr << nbasis << " Non conforming arrays in VectorToMultimap: confs (estimated) num strings is " << confs.size()/nbasis;
    cerr << " isym size is " << isym.size() << " based on nbasis = " << nbasis << endl;
    GA::Terminate();
  }
  
// We need to chop up confs into nbasis chunks of shorts then convert to a string and insert multimap
// Step one skip every nbasis: build an inner vector of length basis inside.
// Easier to use [] operator rather than iterators
// We only need 'sym' every nbasis chunks
  
  int index_sym=-1;
  for (int i = 0; i < confs.size(); i += nbasis ) {
    int sym = isym[++index_sym];
    vector<short> vword; //Assemble a vector of nbasis occupations
    for (int j=0; j< nbasis; ++j) {
      vword.push_back( confs[i+j] );
    }
    string configuration; // Combine all into a string
    convertShortArrayToString( vword, configuration );
    configs.insert(pair<const int,string>(sym,configuration));
  }
#if 0
  cout << "Conversion to Multimap: Size of newly constructed configs multimap is " << configs.size() << endl;
#endif
  return;
}
//Conversion of presumably MPI sent vector of data back into a local multimap array
/* An experimental method to allow scramblinhg of the configutations
   The though is that bunching up by jsym is messing up load-balance with large values of 
   nxtval_chunk

   So instead of using fetchConfig and multimaps, use fetchConfigAlt and Vector-style maps
   Similar code is used for direct-determinant code
*/
void PsociConfigs::convertVectorToVectormap( vector<int> & isym, vector<short> & confs )
{
  cout << "CONVERT TO MAAP" << endl;
  if ( GA::nodeid() == 0 ) cout << "At convertVectorToVectormap method " << endl;
  int nbasis = n_basisftns;
  if ( nbasis == 0 ) { 
    cerr << "convertVectorToVectormap: nbasis is zero cannot recover " << endl;
    GA::Terminate();
  }
  if ( isym.size() == 0 || confs.size() == 0) {
    cout << "Warning: convertVectorToVectorrmap: expected non-empty vector objects. Need to abort" << endl;
    GA::Terminate();
  }
  //Clear out the existing configs vector;
  vec_configs.clear();
  
  // Are confs and isym the same length?
  if ( confs.size()/nbasis != isym.size() ) {
    cerr << nbasis << " convertVectorToVectormap: Non conforming confs and sym arrays: confs (estimated) is of size " << confs.size()/nbasis;
    cerr << " and sym is of size " << isym.size() << " Cannot recover " << endl;
    return;
  }
  if ( confs.size()/nbasis != isym.size() ) {
    cerr << nbasis << " Non conforming arrays in VectorToVectormap: confs (estimated) num strings is " << confs.size()/nbasis;
    cerr << " isym size is " << isym.size() << " based on nbasis = " << nbasis << endl;
    GA::Terminate();
  }
  
// We need to chop up confs into nbasis chunks of shorts then convert to a string and insert multimap
// We only need 'sym' every nbasis chunks
  
  int index_sym=-1;
  for (int i = 0; i < confs.size(); i += nbasis ) {
    int sym = isym[++index_sym];
    vector<short> vword; //Assemble a vector of nbasis occupations
    for (int j=0; j< nbasis; ++j) {
      vword.push_back( confs[i+j] );
    }
    string configuration; // Combine all into a string
    convertShortArrayToString( vword, configuration );
    vec_configs.push_back( pair<const int,string>(sym,configuration) );
  }
#if 0
  cout << "Conversion to Vectormap: Size of newly constructed configs vectormap is " << vec_configs.size() << endl;
#endif

  return;
}

//Note this vector array is the ENTIRE array in a single dimension (nbf * short * numcfgs)
// This may not always work
void PsociConfigs::convertStringToShortArray( string word, vector<short> & vector )
{
  int wordSize = word.size();
  if ( wordSize != n_basisftns ) {
    cerr << "Error in convert to shorts: wrong size " << wordSize << endl;
    cerr << "Aborting " << endl;
    GA::Terminate();
  }
  for (int i=0; i< wordSize; ++i ) {
    vector.push_back( word[i] );
  }
}

//Note this vector array is the ENTIRE array in a single dimension (nbf * short * numcfgs)
void PsociConfigs::convertStringToIntArray( string word, vector<int> & vector )
{
  int wordSize = word.size();
  if ( wordSize != n_basisftns ) {
    cerr << "Error in string to ints: wrong size " << wordSize << endl;
    cerr << "Aborting " << endl;
    GA::Terminate();
  }

// slight modification per akp.
  char element[2];
  for (int i=0; i< wordSize; ++i ) {
    element[0] = word[i];
    element[1] = 0; // guarantees NUL even with optimization
    vector.push_back( atoi(element) );
  }
}

//Note this vector array is the ENTIRE array in a single dimension (nbf * short * numcfgs)
/* This is an alternative: we found weird things happening on Ranger with intel12 
   using the original method
*/
void PsociConfigs::convertStringToIntArray( vector<char> & word, vector<int> & vector )
{
  int wordSize = word.size();
  if ( wordSize != n_basisftns ) {
    cerr << "Error in string to ints: wrong size " << wordSize << endl;
    cerr << "Aborting " << endl;
    GA::Terminate();
  }
  for (int i=0; i< wordSize; ++i ) {
    char element = word[i];
    vector.push_back( atoi(&element) );
  }
}


void PsociConfigs::convertShortArrayToString( vector<short> & vec, string & word )
{
  int wordSize = vec.size();
  if ( wordSize != n_basisftns ) {
    cerr << "Error in convert shorts to string: wrong size " << wordSize << endl;
    cerr << "Aborting " << endl;
    GA::Terminate();
  }
  if ( word.size() > 0 ) {
    cerr << "Warning in convert to string: non-empty string word. Will empty and proceed " << endl;
  }
  word.clear();
  
  vector<short>::iterator it;
  for (it=vec.begin(); it != vec.end(); ++it) {
    //Convert (*it) a short to a char....this should work since the int is 0,1,2 only
    word.push_back( (char)(*it) );
  }
  //cout << "Size of created word is " << word.size() << endl;
}

void PsociConfigs::convertIntArrayToString( vector<int> & vec, string & word )
{
  int wordSize = vec.size();
  if ( wordSize != n_basisftns ) {
    cerr << GA::nodeid() << " Error in convert int to string: wrong size " << wordSize << "n_basisis " <<  n_basisftns<< endl;
    cerr << "Aborting " << endl;
    GA::Terminate();
  }
  if ( word.size() > 0 ) {
    cerr << "Warning in convert to string: non-empty string word. Will empty and proceed " << endl;
  }
  word.clear();
  
  vector<int>::iterator it;
  for (it=vec.begin(); it != vec.end(); ++it) {
    int i = (*it);
    string s;
    stringstream out;
    out << i;
    s = out.str();    
    word.push_back( s[0] );
  }
}

//Can be called locally but nodes should never get out of sync. If one node checks all nodes should consider
//checking.

void PsociConfigs::validateConfigurationSymmetry()
{
  cout << "validateConfigurationSymmetry " << endl;
  cout << " No checking has been implemented at this time " << endl;
}

//Grab the iconfig-th string from the multimap ( convert to C indexing ) 
//This is not so easy becasue of the structure of the maps.... We need another way to do this
//in a distributed environment. We do know, howver, the total number of elements via multimap.size()
//which everyone has.

//Grab all strings s; (start_iconf <= s <= stop_iconf ) INCLUSIVE in some linear order

int PsociConfigs::fetchConfig( int start_iconf, int stop_iconf, vector<pair<int,string> > & selcfgs )
{
  int found=0; // A completely valid response. Real errors trigger a terminate()
  
#ifdef DETAILEDCHECK
  if ( selcfgs.size() != 0 ) {
    cout << "Warning: fetchConfig: expected an empty selcfgs but it was not. Emptying " << endl;
  }
  
  if ( stop_iconf > n_configurations)  {
    cerr << "fetchConfig: An inconsistent specification for iconfig was provided: iconfig = " << stop_iconf;
    cerr << " total number of configurations is currently " << n_configurations << " Aborting: I am " << GA::nodeid()  << endl;
    GA::Terminate();
  }
  
  if ( start_iconf > stop_iconf ) {
    cerr << "fetchConfig: start_iconf >= stop_iconf: " << start_iconf << " " << stop_iconf << " Aborting " << endl;
    GA::error("fetchConfig: start_iconf >= stop_iconf: ", stop_iconf-start_iconf);
  }
#endif

  selcfgs.clear();

  int icount=0;
  multimap<int,string>::iterator it;
  it = configs.begin();
  
  for (it = configs.begin(); it != configs.end(); ++it) {
    ++icount;
    if (icount >= start_iconf && icount <= stop_iconf ) {
      selcfgs.push_back(pair<int,string>( (*it).first, (*it).second ) );
    }
  }

  found = selcfgs.size();
  return( found );
}

/* Doesn't actually get data from GA; that was alrady performed
  process local configs object 
*/
int PsociConfigs::fetchConfigGA(vector<pair<int,string> > & selcfgs, vector<int> & index )
{
  int found=0; // A completely valid response. Real errors trigger a terminate()
  
#ifdef DETAILEDCHECK
  if ( selcfgs.size() != 0 ) {
    cout << "Warning: fetchConfigGA: expected an empty selcfgs but it was not. Emptying " << endl;
  }
#endif

  selcfgs.clear();
  index.clear();

  if ( vec_configs_index.size() != vec_configs.size()  ) GA::error("NEW fetchConfigGA wrong sizes ", vec_configs.size() );

  int size = vec_configs.size();
  for(int i=0; i< size; ++i ) {
   selcfgs.push_back(pair<int,string>( vec_configs[i].first, vec_configs[i].second ) );
   index.push_back( vec_configs_index[i] );
  }

  found = selcfgs.size();
  return( found );
}

int PsociConfigs::fetchConfigAlt( int start_iconf, int stop_iconf, vector<pair<int,string> > & selcfgs )
{
  
#ifdef DETAILEDCHECK
  if ( selcfgs.size() != 0 ) {
    cout << "Warning: fetchConfig: expected an empty selcfgs but it was not. Emptying " << endl;
  }
  
  if ( stop_iconf > n_configurations)  {
    cerr << "fetchConfig: An inconsistent specification for iconfig was provided: iconfig = " << stop_iconf;
    cerr << " total number of configurations is currently " << n_configurations << " Aborting: I am " << GA::nodeid()  << endl;
    GA::Terminate();
  }
  
  if ( start_iconf > stop_iconf ) {
    cerr << "fetchConfig: start_iconf >= stop_iconf: " << start_iconf << " " << stop_iconf << " Aborting " << endl;
    GA::error("fetchConfig: start_iconf >= stop_iconf: ", stop_iconf-start_iconf);
  }
#endif
  
  selcfgs.clear(); // created anew in the calling app
  
  for(int i=start_iconf; i<=stop_iconf;++i) {
    selcfgs.push_back( pair<int,string>( vec_configs[i-1].first, vec_configs[i-1].second ) );
  }
  return( selcfgs.size() );
}

/* Use a primitive shuffle to break up symmetry ordering
   Requires my own rand function etc otherwise different nodes could 
   use different seeds
 
   This doesn't help too much


*/
void PsociConfigs::shuffleConfigurations(int root, vector<int> & isym, vector<short> & confs )
{
  int nbasis = n_basisftns;
  if (GA::nodeid() != root ) return;
  
/* We have a problem. Confs is a contiguous list of shorts with no delimiter. We need to make our own list
 */
  const uint FAKE_SEED = 1;
  srand ( FAKE_SEED ); // This is all we need to guarantee accross nodes the same shuffle
  random_shuffle( isym.begin(), isym.end() );
  
  /* Nw we need to process confs.  How many confs do wew have? ( each of length nbasis )
     Create a new object that is sortable then compact back to confs
  */
  int numconfs = confs.size() / nbasis;
  vector<vector<short> > short_conf; 
  
  for (int i = 0; i < confs.size(); i += nbasis ) {
    vector<short> vword; //Assemble a vector of nbasis occupations
    for (int j=0; j< nbasis; ++j) {
      vword.push_back( confs[i+j] );
    }
    short_conf.push_back( vword ); 
  }
  
  // Now sort
  srand ( FAKE_SEED ); // Use same seed to ensure consistent shuffling
  random_shuffle( short_conf.begin(), short_conf.end() ); //This will shuffle blocks of shorts
  
  // Now push back into configs as one contiguous array
  int index = 0;
  vector<vector<short> >::iterator it;
  vector<short>::iterator sit;
  
  for(it = short_conf.begin(); it != short_conf.end(); ++it ) 
    {
      vector<short> temp = (*it);
      for(sit = temp.begin(); sit != temp.end(); ++sit )
	{
	  confs[index] = (*sit);
	  ++index;
	}
    }
  return;
}

/* Configuation is to be set up as ROWS of length num_basis. Add one extra value to carry the symmetry

*/
void PsociConfigs::createGAconfigs( int num_basis, int num_configs )
{
  const int DATA_TYPE = C_INT;
  const int DATA_NDIM = 2;
  int dims[ 2 ];
  dims[1] = num_configs;
  dims[0] = num_basis+2; // add one for symetry index add another for abs index
  
  const int DISTRIBUTE_EVENLY = -1; //TODO need to optimize this at some point
  
  int chunk[2];
  chunk[0] = dims[0]; //We want distribution by rows (length)
  chunk[1] = DISTRIBUTE_EVENLY;
  
  char local_title[] = "configurations value";
  g_configs = GA::createGA( DATA_TYPE, DATA_NDIM , dims, (char *) local_title, chunk);
  
  char handleGAMessage[] = "PsociPsociConfigs::createGAconfigs";
  g_configs->checkHandle( handleGAMessage);
  g_configs->zero();
  
#if 0
  g_configs->printDistribution();
#endif  
}

void PsociConfigs::destroyGAconfigs()
{
  g_configs->destroy();
}

/* push to GA at location indexGA
   Be sure to APPEND the value of isym to the configuation
   
   Also append the indexGA+1 value, that way subsequent gets can identify the
   absolute indexing of these configurations
   
   LOCAL call BUT intended only for rank=0
*/
void  PsociConfigs::pushConfigsToGA( int indexGA, int nbf, int isym, vector<int> & configs )
{
  if (GA::nodeid() != 0 ) return;
  
  int indexShift = indexGA+1; 
  int n = 1;
  int lo[2];
  lo[1] = indexGA; // used to assist later decision making 
  lo[0] = 0; // 
  int hi[2];
  hi[1] = indexGA;
  hi[0] = nbf-1; // Uhg,these differing index approaches are giving me a headache
  
#if 0
  cout << "PUT ";
  for(int i=0; i< nbf; ++i ) {
    cout << configs[i] ;
  }
  cout << endl;
#endif
  
  g_configs->put( lo, hi, &configs[0], &n);
  
  // now stuff last isym value
  lo[0] = nbf;
  hi[0] = nbf;
  g_configs->put( lo, hi, &isym, &n);
  
  // now stuff last indexGA+1 value
  lo[0] = nbf+1;
  hi[0] = nbf+1;
  g_configs->put( lo, hi, &indexShift, &n);
}


/* Grab all configurations that are LOCAL and push them into a configs
   object

   NOTE at the conclusion, individual cores wil now have (probably) non-zero 
   configs objects, BUT, they will be disjoint. THis is incontrast to the original
   implementation which had configs replicated to all cores

   THe intent being that each node will subsequently convert their set of configs
   to determinants for upload to GA

   convert all distributed calls to vec_configs instrad of straight configs
   LOCAL calls
*/
int PsociConfigs::getLocalConfigsFromGA( )
{
  configs.clear(); //just to be sure
  vec_configs_index.clear();
  vec_configs.clear();
  
  int howmany = 0;
  
  char handleGAMessage[] = "getLocalConfigsFromGA";
  g_configs->checkHandle( handleGAMessage);
  
  int lo[2], hi[2];
  int index;
  int isym;
  
  g_configs->distribution( GA::nodeid(), lo, hi );
  
  if ( lo[0] < 0 || lo[1] < 0 ) {
    cerr << "configuration NO DATA ON ME " << GA::nodeid() << endl;
    return( howmany );
  }
  
  howmany = hi[1] - lo[1] + 1;
  int l_nbf = hi[0] - lo[0] + 1; 
  int n = 1;
  
  int istart = lo[1];
  int istop = hi[1] + 1;
  
  for(int i=istart; i < istop ;  ++i ) {
    vector<int> alllocal( l_nbf );
    
    lo[1]=i; // needs to be an absolute nummber
    hi[1]=i;
    
    g_configs->get( lo, hi, &alllocal[0], &n );

    isym = alllocal[ l_nbf - 2 ]; // we could simply POP
    index = alllocal[ l_nbf - 1 ];
    alllocal.resize( l_nbf - 2 );
    
    string configuration;
    convertIntArrayToString( alllocal, configuration);
    
    vec_configs.push_back( pair<const int,string>(isym, configuration));
    vec_configs_index.push_back( index );
    
    if ( vec_configs_index.size() != vec_configs.size() ) GA::error(" vec_index and vec_config sizes differ" ,configs.size());
  }
  return( howmany );
}


/* return the size of the local configs object */
int PsociConfigs::localConfigsSize()
{
  return ( configs.size() ) ;
}

int PsociConfigs::localVecConfigsSize()
{
  return ( vec_configs.size() ) ;
}





