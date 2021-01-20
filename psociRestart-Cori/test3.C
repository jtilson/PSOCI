/**************************************************************************************
 * Copyright (c) 2011 RENCI.
 * All rights reserved. This program and the accompanying materials
 * MAY BE available under the terms of the RENCI Open Source License
 * UNC at Chapel Hill which accompanies this distribution, and is available at
 * http://www.renci.org/resources/open-source-software-license


**************************************************************************************/
/**
 *
 */
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include <cmath>

int main(int argc, char **argv)
{
   int numit = 128;
   double test;
   int vecsize = 1024 * 1024;

   for(int i=0; i< numit; ++i ) {
      vector<int> testint(vecsize);
      for(int j=0; j< 10; ++j ) {
         vector<int> testintj( 128 );
         test=j*1.0;
      }
   }

}

