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
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2009.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )
*/

/** Partanalyzer
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2009.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )

Compile Options:

g++ -o partanalyzer partanalyzer.cc
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


#ifdef _PARTANALYZER_MAIN
const char* VERSION="alpha 0.5.0.3";

bool 	VERBOSE=false;

bool    QUIET=false;
#endif

//void readListInputFiles(ifstream is, vector<Charr> infilenames){
inline void readListInputFiles(char* argv0, vector<Charr>& infilenames){
        ifstream is(argv0);
        if(!is){
                cout<<"ERROR: Cannot open file "<<argv0<<endl;
                exit(1);
        }
        string fn;
        Charr f;
        if(!QUIET)cout<<"#Getting list of partitions from "<<argv0<<endl;
        while(is>>fn){
                if(fn.compare(0,1,"#")==0){ // It's a comment line...
                        getline(is,fn);
                        continue;               // ignore the whole line and get to the next one.
                }
                f.car = new char[fn.size()+1];
                strcpy(f.car,fn.c_str());
                infilenames.push_back(f);
        }
        if(!QUIET)cout<<"#Partitions list file contains "<<infilenames.size()<<" entries.\n#First entry seen "<<infilenames[0].car<<"\n#Last ("<<infilenames.size()<<") entry seen "<<infilenames[infilenames.size()-1].car<<endl;
}

#endif //END _PARTANALYZER_MAIN_HEADER
