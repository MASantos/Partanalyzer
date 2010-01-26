/*
*    This file is part of Partanalyzer.
*
*    Partanalyzer is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    Partanalyzer is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with Partanalyzer.  If not, see <http://www.gnu.org/licenses/>.
*
*/
/** Partanalyzer
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2011.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
*/

/** Partanalyzer
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2011.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )

Compile Options:

g++ -o partanalyzer partanalyzer.cc
*/

#ifndef _PARTANALYZER_INCLUDES_HEADER
#define _PARTANALYZER_INCLUDES_HEADER 1

#include<iostream>
#include<fstream>
#include<vector>
#include<string>
/**strcmp works well under g++ (GCC) 3.3.3 (SuSE Linux) or g++ (GCC) 4.2.4 (Ubuntu 4.2.4-1ubuntu3 -Gutsy)
but seems to fail to compile under Ubuntu 8.10 Intrepid, 
requiring an explicit inclusion of the C header. An apparent wrapper around it is
the C++ cstring header
*/
#include<cstring>
/**The same happens with
#include<algorithm> //for sort
#include<stdlib.h> //atoi, etc.
#include<string.h> //for strcmp, etc.
*/
#include<algorithm> //for sort
#include<stdlib.h> //atoi, etc.
//
#include<sstream>
#include<map>
#include<set>
#include<list>
#include<math.h>
#include<limits>
#include<time.h>
///Creating , reading and removing directories
#include<dirent.h> ///It's become more or less a standard, although not included in ANSI C(++)
#include<sys/stat.h> ///More random location: others refer to it as living under dir.h or unistd.h ... Looks like rmdir is homeless.
#include<errno.h>
//
/* TESTING CUSTOM GRAPH CLASS 
*/
//#include "graph.h"

using namespace std;

#endif
