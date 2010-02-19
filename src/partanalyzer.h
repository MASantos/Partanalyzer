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
/**@mainpage Partanalyzer
@author Miguel A. Santos, msantos@wodaklab.org
@date 2008-2011
@note Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2011.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )

@detail partanalyzer aims at being a general program for analyzing (sets of) partitions.
 Here a partition is defined as in set theory of mathematics (see
 http://en.wikipedia.org/wiki/Partition_of_a_set). It also allows to
 edit (rudimentarily), as well as generate, partitions.                

*/

#ifndef _PARTANALYZER_MAIN_HEADER
#define _PARTANALYZER_MAIN_HEADER 1

#include "partanalyzer_includes.h"
#include "partanalyzer_definitions.h"
#include "partanalyzer_basic_operations.h"
#include "partanalyzer_help.h"
#include "Roulette.h"
#include "BellNumber.h"
#include "Sequence.h"
#include "MultipleSeqAlign.h"
#include "sNeighborhood.h"
#include "Partition.h" 
#include "PartitionStats.h"
#include "MatrixOfValues.h"
#include "Ccop.h"
#include "RobustDivisiveClustering.h"


#ifdef _PARTANALYZER_MAIN
const char* VERSION="alpha 1.0.";

bool 	DEBUG=false;

bool 	VERBOSE=false;

bool    QUIET=false;

bool	INFO=false;

bool	FUZZYPARTITION=false;
#endif

template<class T> void readListFromFile(char* argv0, T& container){
        if(!QUIET)cout<<"#Reading items list from "<<argv0<<endl;
        ifstream is(argv0);
        if(!is){
                cout<<"ERROR: Cannot open file "<<argv0<<endl;
                exit(1);
        }
        string fn;
        while(is>>fn){
                if(fn.compare(0,1,"#")==0){ // It's a comment line...
                        getline(is,fn);
                        continue;               // ignore the whole line and get to the next one.
                }
		container.push_back(fn);
        }
        if(!QUIET)cout<<"#List file contains "<<container.size()<<" entries.\n#First entry seen "<<container[0]<<"\n#Last ("<<container.size()<<") entry seen "<<container[container.size()-1]<<endl;
}

#endif //END _PARTANALYZER_MAIN_HEADER
