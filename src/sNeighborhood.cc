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

#ifndef _CLASS_SNEIGHBORHOOD
#define _CLASS_SNEIGHBORHOOD 1

#include "sNeighborhood.h"

sNeighborhood::sNeighborhood(){
	point="";
	eps=0.;
	multimode=0;
}
sNeighborhood::sNeighborhood(string p, svect* cl){
	point=p;
        for(cluster::iterator it=cl->begin();it!=cl->end();it++)
                neighbors.insert(*it);
	eps=1.;
	//neighborhoods.insert(pair<sset,double> (neighbors,eps) );
	neighborhoods[neighbors]=eps;
	multimode=0;
}
void sNeighborhood::setMostFrequentNeighborhood(){
	for(map<sset,double>::iterator s=neighborhoods.begin();s!=neighborhoods.end();s++){
		if(s->second==eps)multimode++;
		if(s->second>eps){
			neighbors=s->first;
			eps=s->second;
		}
	}
#ifdef DEBUG
	cout<<"#DEBUG: item "<<point<<" most freq. cluster ("<<eps<<") ";
	for(sset::iterator it=neighbors.begin();it!=neighbors.end();it++)
		cout<<*it<<" ";
	cout<<" Multimode="<<multimode<<endl;
#endif
}

void sNeighborhood::gNeighborsList(){
	for(map<sset,double>::iterator sit=neighborhoods.begin();sit!=neighborhoods.end();sit++){
		for(sset::iterator it=sit->first.begin();it!=sit->first.end();it++)
			cout<<*it<<" ";
		cout<<" ("<<sit->second<<")/ ";
	}
	cout<<" Multimode="<<multimode<<endl;
}

sNeighborhood& sNeighborhood::operator+=(sNeighborhood& snbh){
	if(point.compare("")==0)point=snbh.point;
	if(!point.compare(snbh.point)==0)return *this;
	for(map<sset,double>::iterator it=snbh.neighborhoods.begin();it!=snbh.neighborhoods.end();it++)
		neighborhoods[it->first]+=it->second;
	setMostFrequentNeighborhood();
	return *this;
}

#endif //END _CLASS_SNEIGHBORHOOD
