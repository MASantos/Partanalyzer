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

#ifndef _CLASS_BELLNUMBER_H
#define _CLASS_BELLNUMBER_H 1


///How many partitions exist for a set X? Bell number.
class BellNumber
{
	double _Bn;
	long int _n;
	double Stirling_comb(long int n, long int k);
public:
	BellNumber();
	//double term(long int n);
	double& operator[] (long int n);
	long int n(){ return _n;}
};

#endif //END _CLASS_BELLNUMBER_H
