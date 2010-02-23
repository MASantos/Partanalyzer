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

/** Partanalyzer
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2011.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )

Compile Options:

g++ -o partanalyzer partanalyzer.cc
*/

#ifndef _PARTANALYZER_BASICOPERATIONS_HEADER
#define _PARTANALYZER_BASICOPERATIONS_HEADER 1

#include <vector>
#include <string>
#include <iostream>
#include "partanalyzer_definitions.h"


///Stream out a vector of strings 
inline ostream& operator<<(ostream& os, const svect& cl){
	os<<"{";
	if(cl.empty()){
		os<<"}";
		return os;
	}
	svect::const_iterator it;
	svect::const_iterator clast=cl.end();
	clast--;
	for(it=cl.begin();it!=clast;it++)
		os<<*it<<",";
	os<<*it<<"}";
	return os;
}
///Stream out a set of strings
inline ostream& operator<<(ostream& os, const sset& cl){
	os<<"{";
	if(cl.empty()){
		os<<"}";
		return os;
	}
	sset::const_iterator it;
	sset::const_iterator clast=cl.end();
	clast--;
	for(it=cl.begin();it!=clast;it++)
		os<<*it<<",";
	os<<*it<<"}";
	return os;
}
///Stream out a set of string vectors
inline ostream& operator<<(ostream& os, const set<svect>& cl){
	os<<"{";
	if(cl.empty()){
		os<<"}";
		return os;
	}
	set< svect>::const_iterator it;
	set< svect>::const_iterator clast=cl.end();
	clast--;
	for(it=cl.begin();it!=clast;it++)
		os<<*it<<",";
	os<<*it<<"}";
	return os;
}
///Stream out a vector of string vectors
inline ostream& operator<<(ostream& os, const vector<svect>& cl){
	os<<"{";
	if(cl.empty()){
		os<<"}";
		return os;
	}
	vector< svect>::const_iterator it;
	vector< svect>::const_iterator clast=cl.end();
	clast--;
	for(it=cl.begin();it!=clast;it++)
		os<<*it<<",";
	os<<*it<<"}";
	return os;
}
/*  
*/
/**Algebra of string vectors
a<=b	a is subset of b
*/
inline bool operator<=(const svect& cla, const svect& clb){
	if(DEBUG)cout<<cla<<"<="<<clb<<"?"; 
	long int itfound=0;
	if(cla.size()>clb.size()){
		if(DEBUG)cout<<"NO"<<endl;
		return false;
	}
	if(DEBUG)cout<<"checking..."; 
	for(svect::const_iterator ita=cla.begin();ita!=cla.end();ita++)
		for(svect::const_iterator itb=clb.begin();itb!=clb.end();itb++){
			//cout<<"CL<= : "<<*ita<<"=="<<*itb<<"?"<<endl; 
			if((*ita).compare(*itb)==0)itfound++;
		}
	if(itfound==cla.size()){
		if(DEBUG)cout<<"YES"<<endl;
		return true;
	}
	if(DEBUG)cout<<"NO"<<endl;
	return false;
}
///a<b	a is a strict subset of b
inline bool operator<(const svect& cla, const svect& clb){
	if(DEBUG)cout<<cla<<"<"<<clb<<"?"<<endl; 
	if(cla.size()<clb.size()&& cla<=clb)return true;
	return false;
}
inline bool operator==(svect& cla, svect& clb);

///svect class compare for defining ordered containers (e.g. sets) of svect with the appropriate order
struct svect_compare{
        bool operator() (const svect& cla, const svect& clb) const {
                if(cla.size()==clb.size() && cla==clb) return false;
                return (cla.size()<=clb.size());
                }
};

/*
///Stream out a set of string vectors
inline ostream& operator<<(ostream& os, set< svect>& U){
	os<<"{";
	if(U.empty()){
		os<<"}";
		return os;
	}
	set< svect>::iterator it;
	set< svect>::iterator Ulast=U.end();
	Ulast--;
	for(it=U.begin();it!=Ulast;it++)
		os<<*it<<",";
	os<<*it<<"}";
	return os;
}
*/
///Stream out a set< ,svect_compare> of string vectors
inline ostream& operator<<(ostream& os, const set< svect,svect_compare>& U){
	os<<"{";
	if(U.empty()){
		os<<"}";
		return os;
	}
	set< svect,svect_compare>::const_iterator it;
	set< svect,svect_compare>::const_iterator Ulast=U.end();
	Ulast--;
	for(it=U.begin();it!=Ulast;it++)
		os<<*it<<",";
	os<<*it<<"}";
	return os;
}
///b>=a	a is a subset of b
inline bool operator>=(svect& cla, svect& clb){
	if(clb<=cla)return true;
	return false;
}
///b>a	a is a strict subset of b
inline bool operator>(svect& cla, svect& clb){
	if(cla.size()>clb.size()&& cla>=clb)return true;
	return false;
}
///a!=b	a and b are different
inline bool operator!=(svect& cla, svect& clb){
	return !(cla==clb);
}
///a==b	a and b are different
inline bool operator==(svect& cla, svect& clb){
	long int itfound=0;
	if(cla.size()!=clb.size())return false;
	for(svect::iterator ita=cla.begin();ita!=cla.end();ita++)
		for(svect::iterator itb=clb.begin();itb!=clb.end();itb++)
			if((*ita).compare(*itb)==0)itfound++;
	if(itfound==cla.size()&&cla.size()==clb.size())return true;
	return false;
}
///Intersection of a and b : a*b
inline svect operator*(svect& cla, svect& clb){ 
	svect iab;
	for(svect::iterator ita=cla.begin();ita!=cla.end();ita++)
		for(svect::iterator itb=clb.begin();itb!=clb.end();itb++)
			if((*ita).compare(*itb)==0)iab.push_back(*ita);
	return iab;
}

///Union of clusters of a and b : a+b
inline svect operator+(svect& cla, svect& clb){ 
	svect uab;
	sset suab;
	for(svect::iterator ita=cla.begin();ita!=cla.end();ita++)
		suab.insert(*ita);
	for(svect::iterator itb=clb.begin();itb!=clb.end();itb++)
		suab.insert(*itb);
	for(sset::iterator it=suab.begin();it!=suab.end();it++)
		uab.push_back(*it);
	return uab;
}
///Increment a with b : a=a+b
inline svect& operator+=(svect& cla, svect& clb){
	return (cla=(cla+clb));
}

///Addition of graphs (matrices). Assumes undirected graph and that elements sorted in the same order.
inline graph operator+(graph& ga, graph& gb){
	graph sg;
	for(graph::iterator ita=ga.begin();ita!=ga.end();ita++){
		string a = (ita->first).first ;
		string b = (ita->first).second ;
		strpair sp (a,b);
		edge v=ga[sp];
		if(gb.find(sp)==gb.end())
			sp = make_pair (b,a);
		if(gb.find(sp)!=gb.end()){
			/*cout<<"ERROR: graph operator+ : unmatched pair ("<<a<<","<<b<<") in graph gb:"<<endl;
			exit(1);
			*/
			v+=gb[sp];
		}
		//sg[sp]=ga[sp]+gb[sp];
		sg[sp]=v;
	}
	return sg;
}
///Add gb to ga
inline graph& operator+=(graph& ga, graph& gb){
	return ga=ga+gb;
}
///Graph (matrix) Multiplication by a scalar (double)
inline graph operator*(double& z, graph& ga){
	graph zg;
	for(graph::iterator ita=ga.begin();ita!=ga.end();ita++)
		zg[ita->first]=z*(ita->second);
	return zg;
}
///Simmetric Graph (matrix) Multiplication by a scalar (double)
inline graph operator*(graph& ga, double& z){
	return (z*ga);
}

///Stream out a graph as a pgm picture file 
inline ostream& operator<<(ostream& os, graph& g){
	/*sset its;
	for(graph::iterator i=g.begin();i!=g.end();i++){
		its.insert((i->first).first);
		its.insert((i->first).second);
	}*/
	svect its;
	string fi;
	sset sit;
	graph::iterator i;
	int n=sit.size();
	for( i=g.begin();i!=g.end();i++){
		fi=(i->first).first;
		sit.insert(fi);
		if(sit.size()>n)its.push_back(fi);
		n=sit.size();
		fi=(i->first).second;
		sit.insert(fi);
		if(sit.size()>n)its.push_back(fi);
		n=sit.size();
	}
	double mxv=1.0;
	double mnv=0.0;
	int grey;
	strpair sp;
	os<<"P2"<<endl;
	os<<its.size()<<" "<<its.size()<<endl;
	os<<_PGM_P2_GRAYSCALE_<<endl;
	//for(sset::iterator i=its.begin();i!=its.end();i++){
		//for(sset::iterator j=its.begin();j!=its.end();j++){
	for(vector<string >::iterator i=its.begin();i!=its.end();i++){
		for(vector<string >::iterator j=its.begin();j!=its.end();j++){
			sp=make_pair (*i,*j);
			if(g.find(sp)==g.end()) {
				sp=make_pair (*j,*i);
				if(g.find(sp)==g.end()){
					
					/*
					cout<<"\nERROR: operator<< : graph unknown pair ("<<*i<<","<<*j<<")"<<endl;	
					exit(1);
					*/
				}
			}
			grey=int( (1.0-g[sp])*_PGM_P2_GRAYSCALE_+0.5 );
			if((*i).compare(*j)==0)grey=0;
			cout<<grey<<" ";
		}
		cout<<endl;
	}
	return os;
}

#endif //END _PARTANALYZER_BASICOPERATIONS_HEADER
