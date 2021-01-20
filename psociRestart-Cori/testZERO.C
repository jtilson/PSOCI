/**************************************************************************************

* Copyright (c) 2010 RENCI.
* All rights reserved. This program and the accompanying materials
* are made available under the terms of the RENCI Open Source License.
* UNC at Chapel Hill which accompanies this distribution, and is available at
* http://www.renci.org/resources/open-source-software-license

 Classes: ConstructList

 Description: 

 Author:

 History:


**************************************************************************************/

#include <cstdlib>
#include <string.h>
#include <iostream>
#include <vector>
#include <limits>

#include <sys/sysinfo.h>

using namespace std;

int main(int argc, char **argv) 
{
    cout << "Size of INT is " << sizeof(int) << endl;
    cout << "Size of LONG is " << sizeof(long) << endl;
    cout << "Size of INT* is " << sizeof(int*) << endl;

    cout << "Epsilon is " << std::numeric_limits<double>::epsilon() << endl;


    double a =10.0;
    double b = 10.0+1.0*std::numeric_limits<double>::epsilon();

    if ( a==b ) cout << "equal " << a<<endl;
    cout << " a is " << a <<endl;
    cout << " b is " << b << endl;



}


