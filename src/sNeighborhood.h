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
#include "partanalyzer_definitions.h"

#ifndef _CLASS_SNEIGHBORHOOD_H
#define _CLASS_SNEIGHBORHOOD_H 1


///It defines the neighborhood of a given point, or ball of radius eps around an element 'x'.
///**
//Allows to deal with (sets of) neighborhood(s) of a given point. In particular, it allows to obtain its most frequent
//neighborhood: given a list of sets containing a point p, what is the one that repeats more often.
//If two or more of them have the same frequency, it chooses the first one seen.
//
//It allows the its construction iteratively given other previously instantiated sNeighborhood's.
//
//Its main members point, neighbors and eps, always contain the most frequent neighborhood, with frequency eps, around 'point' (or the largest
//ball around 'point').
//*/
class sNeighborhood
{
public:
	sNeighborhood();
	///Basic constructor: given an element and a set, take the latter as the ball (interval) around that element.
	/**An alternative way of setting up a sNeighborhood is iteratively via the operator += and given other previously
	instantiated sNeighborhood's
	*/
	sNeighborhood(string p, svect* cl);
	//sNeighborhood(neighborhood nb);
	string point;
	///The set of elements of a given ball (interval) around point 'point'.
	sset neighbors;
	///Radius of a given ball (interval) around point 'point'. It's also used as frequency of a given set of observed neighbors.
	double eps;
	///This is strictly speaking, a list of balls (intervals) 'sset' of radii 'double', assumed around point 'point'.
	map<sset,double> neighborhoods;
	///Checks if there are >1 top ranking neighborhoods
	long int multimode;		
	///Determines which neighborhood is the most frequent one.
	/**Runs through each neighborhood in neighborhoods and checks two things:
		if its radius == eps , i.e., to present one, updates multimode++ 
		if radius > eps, selects that neighborhood as the current one, updating also the value of eps to radius.
	*/
	void setMostFrequentNeighborhood();
	//sset::iterator gbegin(){ return (neighborhoods.begin()->first).begin();}
	//sset::iterator gend(){ return (neighborhoods.rbegin()->first).end();}
	sset::iterator lbegin(){ return (neighbors.begin());}
	sset::iterator lend(){ return (neighbors.end());}
	///Prints the list of neighborhoods and correponding radii
	void gNeighborsList();
	///Update list of neighborhoods: adds lists of neighborhoods in snbh to the current one.
	/**This is a local operation, i.e., both sNeighborhoods must be defined around the same point 'point', otherwise it leaves the current one as it is.
	In case the current has no point defined (weird, but lets consider it) it assumes the same one as in snbh.
	Then, it runs over all sets of neighbors in snbh and adds each of them to the current list, neighborhoods. If the one found is equal to
	one already present, both radii are added up; otherwise, the radius is that of the new found set.
	Finally, it sets the most frequent neighborhood as the current one, i.e., that given by, (string) 'point', (double) 'eps' and (sset) 'neighbors'.
	This acts, the facto, as an alterantive constructor of sNeighborhood.
	*/
	sNeighborhood& operator+=(sNeighborhood& snbh);
};

#endif //END _CLASS_SNEIGHBORHOOD_H
