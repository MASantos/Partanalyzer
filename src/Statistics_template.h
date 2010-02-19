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

#ifndef _CLASS_STATISTICS_H
#define _CLASS_STATISTICS_H 1

#include <math.h>
#include <limits>

#define STATISTICS_DEFAULT_NUMBER_VARIABLES 8
///Stores sample values and allows to retrieve them back
/**Stores sample values
*/
template <class T=double>
class Sampling
{
	/**_X contains the _N basic statistics in the following order
	0 - Median
	1 - Mean
	2 - Standard deviation
	3 - Mean error
	4 - Variance
	5 - Mininum value
	6 - Maximum value
	7 - Sample size
	*/
	vector<T> _X;
	int _N;
public:
	Sampling();
	void reset();	
	T median(T v);
	T median();
	T mean(T v);
	T mean();
	T standard_deviation(T v);
	T standard_deviation();
	T mean_error(T v);
	T mean_error();
	T variance(T v);
	T variance();
	T minimum_value(T v);
	T minimum_value();
	T maximum_value(T v);
	T maximum_value();
	T sample_size(T v);
	T sample_size();
	//friend ostream& operator<<(ostream& os, Sampling<T>& S);
	//ostream& operator<<(ostream& os);
};


template<class T>
Sampling<T>::Sampling(){
	reset();
}

template<class T>
void Sampling<T>::reset(){
	_N=STATISTICS_DEFAULT_NUMBER_VARIABLES;
	_X.assign(_N, (T) 0 );
}

template<class T> T Sampling<T>::median(T v){ T ov=_X[0] ; _X[0]=v ; return ov;}
template<class T> T Sampling<T>::median(){ return _X[0];}
template<class T> T Sampling<T>::mean(T v){ T ov=_X[1] ; _X[1]=v ; return ov;}
template<class T> T Sampling<T>::mean(){ return _X[1];}
template<class T> T Sampling<T>::standard_deviation(T v){ T ov=_X[2] ; _X[2]=v ; return ov;}
template<class T> T Sampling<T>::standard_deviation(){ return _X[2];}
template<class T> T Sampling<T>::mean_error(T v){ T ov=_X[3] ; _X[3]=v ; return ov;}
template<class T> T Sampling<T>::mean_error(){ return _X[3];}
template<class T> T Sampling<T>::variance(T v){ T ov=_X[4] ; _X[4]=v ; return ov;}
template<class T> T Sampling<T>::variance(){ return _X[4];}
template<class T> T Sampling<T>::minimum_value(T v){ T ov=_X[5] ; _X[5]=v ; return ov;}
template<class T> T Sampling<T>::minimum_value(){ return _X[5];}
template<class T> T Sampling<T>::maximum_value(T v){ T ov=_X[6] ; _X[6]=v ; return ov;}
template<class T> T Sampling<T>::maximum_value(){ return _X[6];}
template<class T> T Sampling<T>::sample_size(T v){ T ov=_X[7] ; _X[7]=v ; return ov;}
template<class T> T Sampling<T>::sample_size(){ return _X[7];}

//ostream& operator<<(ostream& os, Sampling<T>& S){
/*
template<class T>
ostream& Sampling<T>::operator<<(ostream& os){
	for(typename vector<T>::iterator x=_X.begin();x!=_X.end();x++){
		os<<*it<<"\t";
	}
	return os;
}
*/
#endif //END _CLASS_STATISTICS_H
