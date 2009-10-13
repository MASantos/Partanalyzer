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

#ifndef _CLASS_MATRIXOFVALUES
#define _CLASS_MATRIXOFVALUES 1

#include "MatrixOfValues.h"

MatrixOfValues::MatrixOfValues(){
	_mxofvf=NULL;
	_items.clear();
	_pairs.clear();
	_graph.clear();
	_mx.clear();
	_nedges=0;
	_nitems=0;
	_Tweight=-696969;
}

MatrixOfValues::MatrixOfValues(char* file){
	_mxofvf=file;
	readMxValues();
}

MatrixOfValues::MatrixOfValues(graph Graph){
	_graph=Graph;
	reset();
}

MatrixOfValues::MatrixOfValues(const MatrixOfValues& MOV){
	_graph=MOV.getGraph();
	reset();
	_mxofvf=MOV.FileName();
	//_mxofvf=NULL;
}

/*
MatrixOfValues MatrixOfValues::operator=(MatrixOfValues& MOV){
	_graph=MOV.getGraph();
	reset();
	_mxofvf=MOV.FileName();
}
*/

void MatrixOfValues::reset(){
	_nedges=_graph.size();
	sset items;
	long int oitemssize=0;
	_Tweight=0;
	_pairs.clear();
	_items.clear();
	_mx.clear();
	_nedges=0;
	for(graph::iterator g=_graph.begin();g!=_graph.end();g++){
		_nedges++;
		edge e=g->second; 
		_Tweight+=e;
		_mx.push_back(e);
		node pa=(g->first).first;
		node pb=(g->first).second;
		items.insert(pa);
		if(items.size()>oitemssize){
			_items.push_back(pa);
			oitemssize=items.size();
		}
		items.insert(pb);
		if(items.size()>oitemssize){
			_items.push_back(pb);
			oitemssize=items.size();
		}
		_pairs.insert(pair<string,string> (pa,pb) );
	}
	_nitems=items.size();
	if(_nitems!=_items.size()){
		cout<<"ERROR: MatrixOfValues::reset() : _nitems!=_items.size()"<<endl;
		exit(1);
	}
}

MatrixOfValues MatrixOfValues::pruneEdgesBelow(float edgethreshold, bool terse){
	long int nprunedges=0;
	return pruneEdges(edgethreshold, nprunedges, true, terse);
}
MatrixOfValues MatrixOfValues::pruneEdgesBelow(float edgethreshold, long int& nprunedges, bool terse){
	return pruneEdges(edgethreshold, nprunedges, true, terse);
}
MatrixOfValues MatrixOfValues::pruneEdgesAbove(float edgethreshold, bool terse){
	long int nprunedges=0;
	return pruneEdges(edgethreshold, nprunedges, false, terse);
}
MatrixOfValues MatrixOfValues::pruneEdgesAbove(float edgethreshold, long int& nprunedges, bool terse){
	return pruneEdges(edgethreshold, nprunedges, false, terse);
}
MatrixOfValues MatrixOfValues::pruneEdges(float edgethreshold, long int& nedgesdel, bool below, bool terse){
	if(terse&&!QUIET){
		string s=below?"below":"above";
		cout<<"#Pruning edges "<<s<<" "<<edgethreshold<<endl;
	}
	long int ndeledges=0;
	///If we erase edges we'll lose the singletons. Thus, instead just label those edges as NaN
	graph prunedGraph(_graph.begin(),_graph.end());
	for(graph::iterator e=_graph.begin();e!=_graph.end();e++){
		if(below){
			if(e->second<edgethreshold){
				//prunedGraph.erase(e->first);
				prunedGraph[e->first]=numeric_limits<edge>::quiet_NaN();
				ndeledges++;
			}
		}
		else if(e->second>=edgethreshold){
				//prunedGraph.erase(e->first);
				prunedGraph[e->first]=numeric_limits<edge>::quiet_NaN();
				ndeledges++;
			}
	}
	if(terse&&!QUIET)cout<<"#Edges deleted= "<<ndeledges<<endl;
	nedgesdel=ndeledges-nedgesdel;
	return MatrixOfValues(prunedGraph);
}

Partition MatrixOfValues::cluster(){
	map<node, sset* > canonicalMapping;
	set<node > coset;
	if(DEBUG)cout<<"#Clustering graph\n#First generate the clusters as ssets (=set<node>)"<<endl;
	///First generate the clusters as ssets (=set<node>)
	for(graph::iterator e=_graph.begin();e!=_graph.end();e++){
		///For each pairs of nodes A and
		node A=(e->first).first;
		///B
		node B=(e->first).second;
		if(DEBUG){
			if(isnan(e->second)){
				cout<<"#Found unconnected nodes {A="<<A<<" , B="<<B<<"}= "<<(e->second)<<endl;
			}else{
				cout<<"#Found edge (A="<<A<<" , B="<<B<<")= "<<(e->second)<<endl;
			}
		}
		if(!(A<coset)){
		/**If A hasn't been already assigned to a class, then
			if B wasn't assigned either
			*/
			if(!(B<coset)){
				///Allocate memory for a new class (A's)
				sset* sptr = new sset ;
				///Equivalence class of A gets its inaugural (& representative) member A
				sptr->insert(A);
				///If there is no actual edge between A and B
				if(isnan(e->second)){
					///Allocate memory for a new class (B's)
					sset* bsptr = new sset ;
					///Equivalence class of B gets its inaugural (& representative) member B
					bsptr->insert(B);
					if(DEBUG)cout<<"#First time seing A and B: Each assigned to different clusters: ("<<sptr<<") "<<*sptr<<" and ("<<bsptr<<") "<<*bsptr<<endl;
					///Store each member's pointer, which point 
					canonicalMapping.insert(pair<node, sset*> (A,sptr));
					///to that same class
					canonicalMapping.insert(pair<node, sset*> (B,bsptr));
				///If there is an actual edge between A and B
				}else{
					///and B is in the same class as A
					sptr->insert(B);
					if(DEBUG)cout<<"#First time seing A and B: Both assigned to the same cluster: ("<<sptr<<") "<<*sptr<<endl;
					///Store each member's pointer, which point 
					canonicalMapping.insert(pair<node, sset*> (A,sptr));
					///to that same class
					canonicalMapping.insert(pair<node, sset*> (B,sptr));
				}
				///Update the list of assigned nodes
				coset.insert(A);
				coset.insert(B);
			}
			///else, if B had already been assigned to a class
			else{
				///If there is no actual edge between A and B
				if(isnan(e->second)){
					///Allocate memory for a new class (A's)
					sset* sptr = new sset ;
					///Equivalence class of A gets its inaugural (& representative) member A
					sptr->insert(A);
					if(DEBUG)cout<<"#First time seing A : assigned to each own cluster: ("<<sptr<<") "<<*sptr<<endl;
					///Store A's own class pointer (pointing to A's class) 
					canonicalMapping.insert(pair<node, sset*> (A,sptr));
				///If there is an actual edge between A and B
				}else{
					///retrieve B's class pointer
					set<node >* bptr=canonicalMapping[B];
					///and insert A in the same class as B 
					bptr->insert(A);
					if(DEBUG)cout<<"#First time seing A : assigned to the same cluster as B: ("<<bptr<<") "<<*bptr<<endl;
					///Store A's class pointer (pointing to A's class) as the same pointer B has
					canonicalMapping.insert(pair<node, sset*> (A,bptr));
				}
				///Update the list of assigned nodes
				coset.insert(A);
			}
		}else{
		/**However, if A was already assigned to a class, then
			if B was not,
			*/
			if(!(B<coset)){
				///If there is no actual edge between A and B
				if(isnan(e->second)){
					///Allocate memory for a new class (B's)
					sset* bsptr = new sset ;
					///Equivalence class of B gets its inaugural (& representative) member B
					bsptr->insert(B);
					if(DEBUG)cout<<"#First time seing B : assigned to each own cluster: ("<<bsptr<<") "<<*bsptr<<endl;
					///Store B's class pointer (pointing to B's class) as the same pointer A has
					canonicalMapping.insert(pair<node, sset*> (B,bsptr));
				///If there is an actual edge between A and B
				}else{
					///retrieve A's class pointer
					set<node >* aptr=canonicalMapping[A];
					///and insert B in the same class as A 
					aptr->insert(B);
					if(DEBUG)cout<<"#First time seing B : assigned to the same cluster as A: ("<<aptr<<") "<<*aptr<<endl;
					///Store B's class pointer (pointing to B's class) as the same pointer A has
					canonicalMapping.insert(pair<node, sset*> (B,aptr));
				}
				///Update the list of assigned clusters
				coset.insert(B);
			}
			/**else, if B had also been already assigned to a class
			we need to merge both clusters (A's and B's) if this is an actual edge, so...
			*/ 
			else{
				///If there is no actual edge between A and B, do nothing, but
				if(isnan(e->second)){
					if(DEBUG)cout<<"#A and B were already assigned, but No edge between A and B: Nothing to do."<<endl;
				///If there is an actual edge between A and B
				}else{
					///retrieve A's class pointer
					set<node >* aptr=canonicalMapping[A];
					///retrieve B's class pointer
					set<node >* bptr=canonicalMapping[B];
					///There shouldn't be any duplicate edge, but just in case, let's not waste resources
					if(aptr==bptr){
						if(DEBUG)cout<<"#A and B belong to the same class. Nothing to do."<<endl;
						continue;
					}else{
						if(DEBUG)cout<<"#A and B were already assigned to, respectively, ("<<aptr<<") "<<*aptr<<" and ("<<bptr<<") "<<*bptr<<endl;
					}
					///merge B into A
					(*aptr)+=(*bptr);
					if(DEBUG)cout<<"#Merging them now into A's cluster : ("<<aptr<<") "<<*aptr<<endl;
					///Update B's class pointer (pointing to B's class), as well as that of each of its relatives, to the same pointer A has
					for(sset::iterator n=bptr->begin();n!=bptr->end();n++) 
						canonicalMapping[*n]=aptr;
					///Delete B's old cluster
					delete bptr;
				}
			}
			
		}
	}
	///Then generate a set of unique sset* pointers. Each point to one of the distinct classes (clusters) defined earlier.
	set<sset* > sclassp;
	for(map<node, sset*>::iterator cm=canonicalMapping.begin();cm!=canonicalMapping.end();cm++){
		sclassp.insert(cm->second);
		if(DEBUG)cout<<"#CanonicalMap: "<<cm->first<<"\t"<<cm->second<<"\t"<<*(cm->second)<<endl;
	}
	if(DEBUG){
		for(set<sset*>::iterator pp=sclassp.begin();pp!=sclassp.end();pp++)
			cout<<"#ClusterFound: "<<*pp<<"\t"<<**pp<<endl;
	}
	int offset=2;
	char* partf=NULL;
	char* tabf=NULL;
	///Finally, generate the (fully-sorted) partition and return it as a value
	return Partition(sclassp,offset,partf,tabf);
}

Sampling MatrixOfValues::edgeDistributionStats(){
	double avg=0;
	double var=0;
	double std=0;
	double stder=0;
	double min=numeric_limits<double>::max();
	double max=numeric_limits<double>::min();
	double v=0;
	double vv=0;
	double x;
	for(graph::iterator g=_graph.begin();g!=_graph.end();g++){
		x=g->second;
		v+=x;
		vv+=(g->second)*(g->second);
		min=(x<min)?x:min;
		max=(x>max)?x:max;
	}
	avg=v/_nedges;
	var=vv/_nedges-pow(avg,2); std=sqrt(var); std*=sqrt(_nedges*1.0/(_nedges-1.0)); //Sample standard deviation
	//
	Sampling stats;
	stder=std/sqrt(1.0*_nedges);
	stats.median(0.0); //0 : Will become the median later
	stats.mean(avg); //1 : Average edge
	stats.standard_deviation(std); //2 : Standard deviation
	stats.mean_error(stder);//3: Mean error
	stats.variance(var); //4 : Variance
	stats.minimum_value(min); //5 : Minimum value
	stats.maximum_value(max); //6 : Maximum value
	stats.sample_size((double)_nedges); //sample size
	return stats;
}
void MatrixOfValues::edgeDistribution(Partition* part){
	sset nodes;
	cout<<"#BeginEdgeDistribution"<<endl;
	cout<<"#Node\t\tAverage\t\tStd.Dev\t\tStd.Err.\tSkew.\t\tMin.\t\tMax.\t\tSampleSize";
	if(!(part==NULL))cout<<"\tClusterSize\tCluster"<<endl;
	else cout<<endl;
	int n=nodes.size();
	int clsize=-1;
	string clname="";
	string tmps;
	for(smap::iterator p=_pairs.begin();p!=_pairs.end();p++){
		nodes.insert(p->first);
		if(n<nodes.size()){
			n=nodes.size();
			if(!(part==NULL)){
				tmps=p->first;
				clsize=part->getClusterSize(p->first);
				clname=part->getClusterName(tmps);
			}
			edgeDistribution(p->first,clsize,clname);
		}
		nodes.insert(p->second);
		if(n<nodes.size()){
			n=nodes.size();
			if(!(part==NULL)){
				tmps=p->second;
				clsize=part->getClusterSize(p->second);
				clname=part->getClusterName(tmps);
			}
			edgeDistribution(p->second,clsize,clname);
		}
	}
	cout<<"#EndEdgeDistribution"<<endl;
}

void MatrixOfValues::edgeDistribution(const string& sa, int clusterSize, string clName){
	double min=9999999999.0;
	double max=-min;
	double avg=0.0;
	double std=0.0;
	double skewness=0.0;
	int samplesize=0;
	//edge v;
	for(graph::iterator g=_graph.begin();g!=_graph.end();g++){
		const edge& v=g->second;
		const string& sb1=(g->first).first;
		const string& sb2=(g->first).second;
		if(sa.compare(sb1)==0||sa.compare(sb2)==0){
			samplesize++;
			avg+=v;
			std+=pow(v,2);
			skewness+=pow(v,3);
			min=v<min?v:min;
			max=v>max?v:max;
		}
	}
	avg/=(1.0*samplesize);
	std/=(1.0*samplesize);
	skewness=skewness/(1.0*samplesize);
	skewness=skewness-3*avg*std+2*pow(avg,3);
	std-=pow(avg,2);
	std=sqrt(std);
	skewness/=pow(std,3);
	printf("%s\t%4.2e\t%4.2e\t%4.2e\t%4.2e\t%4.2e\t%4.2e\t%d",sa.c_str(),avg,std,std/sqrt(1.0*samplesize),skewness,min,max,samplesize);
	if(clusterSize>0)printf("\t%d\t%s\n",clusterSize,clName.c_str());
	else printf("\n");
	//cout<<sa<<"\t"<<avg<<"\t"<<std<<"\t"<<std/sqrt(1.0*samplesize)<<"\t"<<skewness<<"\t"<<min<<"\t"<<max<<"\t"<<samplesize<<endl;
}

strpair MatrixOfValues::cullEdge(string a, string b){
	if(_graph.find(strpair (a,b))!=_graph.end() ) return strpair (a,b);
	if(_graph.find(strpair (b,a))!=_graph.end() ) return strpair (b,a);
	cout<<"ERROR: MatrixOfValues::cullEdge : "<<FileName()<<" : Not found edge between nodes "<<a<<" and "<<b<<endl;
	exit(1);
}

bool MatrixOfValues::existEdge(string a, string b){
	if(_graph.find(strpair (a,b))!=_graph.end() || _graph.find(strpair (b,a))!=_graph.end() )
		return true;
	else
		return false;
}

double MatrixOfValues::v(int i, int j) {
	int k;
	k=((2*_nitems-2-i)*(i+1))/2+j-_nitems;
	//cout<<"k... "<<(_nitems-1-i/2.0)*(i+1)<<"+"<<j<<"-"<<_nitems<<endl;
	//cout<<"mx["<<k<<"]="<<_mx[k]<<endl;
	return _mx[k] ;
}
double MatrixOfValues::v(string a , string b, REDMxVal useRED){ //useRED Defaults to useOrgRED
	/**As we are passing a copy, we can change here a,b without changing the orignal variables.
	If it's an aritifically generated redundant sequences (-RED-#), use same matrix values as original. 
	Drop substring tail starting at -RED
	If useOrgRED, matrix value for redundant sequence is that of the original one: IDEAL DUPLICATE.
	If not, see if useOwnRED: redundant sequences are expected to have their own defined matrix values, e.g., introducing noise on the origingal
	matrix.
	if it is useZeroRED, use default edge value (see below) for the redundant sequences.
	Otherwise, not yet defined. Meanwhile, like useZeroRED.
	*/
	switch (useRED ){
		case useOrgRED:
				a=a.substr(0,a.find("-RED"));
				b=b.substr(0,b.find("-RED"));
				///Should be valid only for the artificially generated duplicates. Care should be taken that it doesn't affect other cases
				if(a.compare(b)==0)return 1.0;
				break;
		case useOwnRED:
		case useZeroRED:
		default:
				break;
	}
	if(VERBOSE)cout<<"#Searching for pair "<<a<<" , "<<b<<endl;
	if(_graph[pair<string,string> (a,b)]) return _graph[pair<string,string> (a,b)];
	if(_graph[pair<string,string> (b,a)]) return _graph[pair<string,string> (b,a)];
	/**Using a map with doubles as values (pair of strings as keys) does not allow in an easy way to distinguish between
	    a key with a value=0 or a key that simply was not previously defined.
	    as we are dealing here with graphs, i.e, those keys represent edges, and the values represent numerical weight 
	    (and given that I haven't found an easy walkaround) for the time being I'll silently ignore and just print the
	    value for the first key. Still I leave the checking above, otherwise we may assign a zero to an otherwise non-zero edge.
	   cout<<"ERROR in reading matrix "<<_mxofvf<<" : Pair not found "<<a<<" , "<<b<<endl; */
	return _graph[pair<string,string> (a,b)];
	exit(1);
}

int MatrixOfValues::_getIndexOfItem(string str){
	int i;
	for(i=0;i<_items.size();i++)
		if( strcmp(str.c_str() , _items[i].c_str())==0 )
			return i;
	cout<<"ERROR: index for "<<str<<" not found. nitems= "<<_items.size()<<" i="<<i<<endl;
	exit(1);
}

void MatrixOfValues::readMxValues(){ ///Later on we'll assume _mx represents a square matrix and we don't care about the diagonal values. 
	string pa,pb,s;	///Thus input matrix element (i,j) is located at index k=(_nitems-1-i/2)*(i+1)+j-_nitems  of vector _mx, where i,j=0,1,2,... and 
	edge v;		///
	long int rep;
        ifstream _is(_mxofvf);
        if(!_is){
                cout<<"Cannot open file "<<string(_mxofvf)<<endl;
                exit(1);
        }
	row mrow;
	_nedges=0;
	_nitems=0;
	_Tweight=0;
	int a,b,c,m;
	if(!QUIET) cout<<"#Reading matrix "<<_mxofvf<<" ... ";
	while(_is>>pa){
		if(pa.substr(0,1).compare("#")==0){
			getline(_is,pa);
			continue;
		}
		_is>>pb ; _is>>v;
#ifdef DEBUGDETAILS
		cout<<"Seen: "<<pa<<" "<<pb<<" "<<v<<endl;
#endif
		if(pa==pb) continue; 
		a=b=c=0;
		for(int i=0;i<_items.size();i++){
			if(strcmp(_items[i].c_str(),pa.c_str())==0)a=1;
			if(strcmp(_items[i].c_str(),pb.c_str())==0)b=1;
			//if(a==1&&b==1)break;
		}
		
		if(a==0){ _items.push_back(pa); _nitems++;}
		if(b==0){ _items.push_back(pb); _nitems++;}
		//if(a*b==1)continue;
		_nedges++;
		_mx.push_back(v);
		_pairs.insert(pair<string,string> (pa,pb) );
		//_graph.insert( pair<pair<string,string>,edge> (pair<string,string> (pa,pb), v) );
		//pair<graph::iterator,bool> ret;
		_graph[pair<string,string> (pa,pb)]=v ;
		_Tweight+=v;
#ifdef DEBUG
		cout<<"#inserting pair: "<<pa<<" "<<pb<<" #pairs="<<_pairs.size()<<endl;
#endif
	}
	if(!QUIET) cout<<"#\n#Finish Reading matrix: #items="<<_nitems<<" #matrixelements="<<_mx.size()<<" #pairs="<<_pairs.size()<<" TWeight="<<_Tweight<<endl;
}

void MatrixOfValues::cull(MatrixOfValues* mx2){
        if(!QUIET) cout<<"#Cull edges specified in list= "<<mx2->_mxofvf<<" from graph mx= "<<_mxofvf<<endl;
        for(smap::iterator it2=mx2->_pairs.begin();it2!=mx2->_pairs.end();it2++){
		strpair sp=cullEdge(it2->first,it2->second);
		cout<<sp.first<<"\t"<<sp.second<<"\t"<<v(sp.first,sp.second)<<endl;
	}
}

void MatrixOfValues::merge(MatrixOfValues* mx2){
        if(_pairs.size()!=mx2->_pairs.size()){
                cout<<program<<" : Error: graphs contain different number of edges : mx1= "<<_pairs.size()<<" mx2= "<<mx2->_pairs.size()<<endl;
                exit(1);
        }
        if(!QUIET) cout<<"#Merging graph mx1="<<_mxofvf<<" + mx2="<<mx2->_mxofvf<<endl;
        for(smap::iterator it=_pairs.begin();it!=_pairs.end();it++)
                cout<<it->first<<"\t"<<it->second<<"\t"<<v(it->first,it->second)<<"\t"<<mx2->v(it->first,it->second)<<endl;
}

void MatrixOfValues::merge(MatrixOfValues* mx2, Partition* pt){
        if(_pairs.size()!=mx2->_pairs.size()){
                cout<<program<<" : Error: graphs contain different number of edges : mx1= "<<_pairs.size()<<" mx2= "<<mx2->_pairs.size()<<endl;
                exit(1);
        }
        if(!QUIET) cout<<"#Merging graph mx1="<<_mxofvf<<" + mx2="<<mx2->_mxofvf<<" and partition "<<pt->FileName()<<endl;
	string clname;
	map<string,int> unassigneditems;
        for(smap::iterator it=_pairs.begin();it!=_pairs.end();it++){
		clname=pt->areWithinSameCluster(it->first,it->second);
		if(strcmp(clname.c_str(),"NAN1")==0)
			unassigneditems.insert( pair<string,int> (it->first,1) );
		else if(strcmp(clname.c_str(),"NAN2")==0)
			unassigneditems.insert( pair<string,int> (it->second,1) );
                //cout<<v(it->first,it->second)<<"\t"<<mx2->v(it->first,it->second)<<"\t"<<clname<<"\t"<<it->first<<"\t"<<it->second<<endl;
                ///Old version of overlapVspearson-toxmgrace script requires an additional dummy column...to be removed soon...
		cout<<v(it->first,it->second)<<"\t"<<mx2->v(it->first,it->second)<<"\t"<<2<<"\t"<<clname<<"\t"<<it->first<<"\t"<<it->second<<endl;
	}
	if(!QUIET) if(unassigneditems.size()>0){
		cout<<"#Warning: items not found in any cluster: "<<unassigneditems.size()<<" ";
		for(map<string,int>::iterator it=unassigneditems.begin();it!=unassigneditems.end();it++)
			cout<<it->first<<" ";
		cout<<endl;
	}
}

void MatrixOfValues::printMatrix(){
	if(!QUIET)cout<<"#items "<<_nitems<<endl;
	if(!QUIET)cout<<"#Number of pairs: "<<_pairs.size()<<endl;
	for(smap::iterator it=_pairs.begin();it!=_pairs.end();it++)
		cout<<it->first<<"\t"<<it->second<<"\t"<<v(it->first,it->second)<<endl;
}

#endif //END _CLASS_MATRIXOFVALUES
