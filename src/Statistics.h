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

#ifndef _CLASS_STATISTICS_H
#define _CLASS_STATISTICS_H 1

#include<vector>
#include<math.h>
#include<limits>
#include<iostream>

using namespace std;

#define STATISTICS_DEFAULT_NUMBER_VARIABLES 8
/**Stores sample values
 * */
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
	vector<double > _X;
	int _N;
public:
	Sampling();
	void reset();	
	double median(double v);
	double median();
	double mean(double v);
	double mean();
	double standard_deviation(double v);
	double standard_deviation();
	double mean_error(double v);
	double mean_error();
	double variance(double v);
	double variance();
	double minimum_value(double v);
	double minimum_value();
	double maximum_value(double v);
	double maximum_value();
	double sample_size(double v);
	double sample_size();
	friend ostream& operator<<(ostream& os, Sampling& S);
};

#endif //END _CLASS_SdoubleAdoubleISdoubleICS_H
