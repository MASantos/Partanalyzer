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

#include "partanalyzer_includes.h"
#include "BellNumber.h"

#ifndef _CLASS_BELLNUMBER
#define _CLASS_BELLNUMBER 1


///How many partitions exist for a set X? Bell number.

BellNumber::BellNumber(){
	_Bn=1;
	_n=0;
}
double BellNumber::Stirling_comb(long int n, long int k){
	double c;
	if(n==k | k==0)return 1.0;
	if(n<k || n<0 || k<0 ){
		cout<<"ERROR: Stirling_comb : n<k | n<0 | k<0"<<endl;
		exit(1);
	}
	c=sqrt((double)n/(double)(k*n-k*k));
	c*=pow((double) n, (double) n);
	c/=(pow((double)k,(double)k)*pow((double)(n-k),(double)(n-k)));
	return c;
}

//double BellNumber::term(long int n=0){
double& BellNumber::operator[](long int n){
	_Bn=0.0;
	BellNumber B;
	if(n==0) _Bn=1;
	for(long int k=0;k<n;k++)
		_Bn+=Stirling_comb(n,k)*B[k];
		//_Bn+=Stirling_comb(n,k)*term(k);
	_n=n;
	return _Bn;
}

#endif //END _BELLNUMBER
