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


#ifndef _CLASS_ROULETTE
#define _CLASS_ROULETTE 1

#include "Roulette.h"

roulette::roulette(){
}
/// n random integers within interval [0,n)
vector< int> roulette::spinWheeli(int n){
	spinWheeli(n,n);
}
/// n random integers within interval [0,L)
vector< int> roulette::spinWheeli(int n, int L){ 
	vector< int > randseq;
	int i=0;
	while(n>i++){
#ifndef DEBUG
		randseq.push_back((int)floor((double)L*rand()*1.0/RAND_MAX));
#else
		int myr=rand();
		double dmyr=(double)L*myr*1.0/RAND_MAX;
		randseq.push_back((int)floor(dmyr));
		cout<<"#DEBUG: roulette : spinWheeli "<<i<<"("<<n<<","<<L<<") : "<<myr<<"("<<RAND_MAX<<","<<dmyr<<")"<<" :  "<<randseq[randseq.size()-1]<<endl;
		getchar();
#endif
	}
	return randseq;
}
///Writes the seed used in the current directory under the file .seed_used
void roulette::writeSeed(){
	const char* dir="./";
	writeSeed(dir);
}
///Writes the seed used in the provided directory under the file .seed_used
void roulette::writeSeed(const char* dir){
	stringstream ss;
	ss<<dir<<"/.seed_used";
	ofstream of(ss.str().c_str());
	if(!of){
		cout<<"#WARNING: coud not open file .seed_used"<<endl;
	}else{
		of<<_seed<<endl;
		of.close();
	}
}

#endif //END _ROULETTE
