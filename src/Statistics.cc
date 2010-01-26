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

#ifndef _CLASS_STATISTICS
#define _CLASS_STATISTICS 1

#include "Statistics.h"

Sampling::Sampling(){
	reset();
}

void Sampling::reset(){
	_N=STATISTICS_DEFAULT_NUMBER_VARIABLES;
	_X.assign(_N, 0.0 );
}

double Sampling::median(double v){ double ov=_X[0] ; _X[0]=v ; return ov;}
double Sampling::median(){ return _X[0];}
double Sampling::mean(double v){ double ov=_X[1] ; _X[1]=v ; return ov;}
double Sampling::mean(){ return _X[1];}
double Sampling::standard_deviation(double v){ double ov=_X[2] ; _X[2]=v ; return ov;}
double Sampling::standard_deviation(){ return _X[2];}
double Sampling::mean_error(double v){ double ov=_X[3] ; _X[3]=v ; return ov;}
double Sampling::mean_error(){ return _X[3];}
double Sampling::variance(double v){ double ov=_X[4] ; _X[4]=v ; return ov;}
double Sampling::variance(){ return _X[4];}
double Sampling::minimum_value(double v){ double ov=_X[5] ; _X[5]=v ; return ov;}
double Sampling::minimum_value(){ return _X[5];}
double Sampling::maximum_value(double v){ double ov=_X[6] ; _X[6]=v ; return ov;}
double Sampling::maximum_value(){ return _X[6];}
double Sampling::sample_size(double v){ double ov=_X[7] ; _X[7]=v ; return ov;}
double Sampling::sample_size(){ return _X[7];}

ostream& operator<<(ostream& os, Sampling& S){
	for(vector<double>::iterator x=S._X.begin();x!=S._X.end();x++){
		os<<*x<<"\t";
	}
	return os;
}

#endif //END _CLASS_STATISTICS
