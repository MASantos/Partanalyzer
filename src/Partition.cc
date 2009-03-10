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

#ifndef _CLASS_PARTITION
#define _CLASS_PARTITION 1

#include "Partition.h"

Partition::Partition(){
}

graph Partition::setAdjacencyMatrix(){
	_Ad.clear();
	if(sitems.empty()) getItems();
	for(sset::iterator ita=sitems.begin();ita!=sitems.end();ita++){
		_Ad[ strpair (*ita,*ita) ]= 1.0;
		if(*ita==*sitems.rbegin()) break; //Otherwise segfault, as j would go past end
		sset::iterator itb=ita; //Can't declare & double initialize within for-loop?
		for(itb++ ;itb!=sitems.end();itb++)
			if( areEquiv(*ita,*itb) ) _Ad[ strpair (*ita,*itb) ]= 1.0;
			else _Ad[ strpair (*ita,*itb) ]= 0.0;
	}
	if(_Ad.empty()){
		cout<<"\nERROR: Partition::setAdjacencyMatrix() : Failed to build Adjacency Matrix : Partition "<<FileName()<<endl;
		exit(1);
	}
	//if(VERBOSE)cout<<"#TEST\n"<<_Ad;
	return _Ad;
}

graph Partition::setAdjacencyMatrix_os(){
	_Ad.clear();
        for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++){
                for(cluster::iterator ita=cl->begin()+_items_offset;ita!=cl->end();ita++){
                        for(cluster::iterator itb=cl->begin()+_items_offset;itb!=cl->end();itb++){
                                _Ad[ strpair (*ita,*itb) ]= 1.0;
                        }
                }
        }
	/*
	for(vector<string >::iterator ita=_vitems.begin();ita!=_vitems.end();ita++){
		_Ad[ strpair (*ita,*ita) ]= 1.0;
		if(*ita==*_vitems.rbegin()) break; //Otherwise segfault, as j would go past end
		vector<string >::iterator itb=ita; //Can't declare & double initialize within for-loop?
		for(itb++ ;itb!=_vitems.end();itb++){
			if( areEquiv(*ita,*itb) ) _Ad[ strpair (*ita,*itb) ]= 1.0;
			else _Ad[ strpair (*ita,*itb) ]= 0.0;
			if(VERBOSE)cout<<*ita<<" "<<*itb<<endl;
		}
	}
	*/
	if(_Ad.empty()){
		cout<<"\nERROR: Partition::setAdjacencyMatrix_os() : Failed to build Adjacency Matrix : Partition "<<FileName()<<endl;
		exit(1);
	}
	if(VERBOSE)cout<<"#TEST\n"<<_Ad;
	return _Ad;
}

/*
Partition::Partition(char* file,partFileFormat iformat=partFmtPART){
	_items_offset=2;
	_partitionf=file;
	_piformat=iformat;
	_readClusters();
}
*/
void Partition::_resetMembers(){
	_nclusters=clusters.size();
	_nitems=_nsingletons=0;
	_largest_cluster=-1;
	_it_largest_cluster=clusters.end();
	sitems.clear();
	ssingletons.clear();
	long int mxcls=0;
	int mxclidx=0;
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++ ,mxclidx++){
		if(cl->size()-_items_offset>mxcls){
			mxcls=cl->size();
			_largest_cluster=mxclidx;
			_it_largest_cluster=cl;
		}
		if(cl->size()-_items_offset==1) {
			_nsingletons++; 		//We might assign its value at the end to ssingletons.size(), but this way allows for checking of duplicate singleton clusters, i.e., non-sound partition.
			ssingletons.insert( (*cl)[_items_offset] );
		}
		for(svect::iterator it=cl->begin()+_items_offset;it!=cl->end();it++ , _nitems++)
			sitems.insert(*it);
	}
	//_nsingletons=ssingletons.size(); ///See comment above
	//setAdjacencyMatrix();
}

Partition::Partition(char* file, partFileFormat iformat, int ofs){//Defaults: iformat=partFmtPART, ofs=2
	_items_offset=ofs;
	_partitionf=file;
	_piformat=iformat;
	_mcltabf=NULL;
	_readClusters();
//	if(_largest_cluster<0||_largest_cluster>=clusters.size())
//		if(!VERBOSE)cout<<"#WARNING: Partition::Partition : out_of_range index _largest_cluster "<<endl;
}

Partition::Partition(smat* clustersl, int ofs, char* partf, char* tabf){//Defaults: partf=NULL, tabf=NULL
	//_partitionf=NULL;
	//_mcltabf=NULL;
	_partitionf=partf;
	_mcltabf=tabf;
	_piformat=partFmtPART;
	clusters=*clustersl;
	_items_offset=ofs;
	_resetMembers();
	/*
	_nclusters=clusters.size();
	_nsingletons=0;
	_nitems=0;
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++){
		if(cl->size()-_items_offset==1){
			_nsingletons++;
			ssingletons.insert(cl->begin()+_items_offset,cl->end());
		}
		_nitems+=cl->size()-_items_offset;
	}
	*/
	if(_nsingletons!=ssingletons.size()){
		cout<<"ERROR: Partition::_readClusters : Not a sound Partition : non mutually disjoint singleton clusters"<<endl;
		cout<<ssingletons<<endl;
		exit(1);
	}
	getItems(); //Update sitems
	if(_nitems!=sitems.size() || ! sitems.size()>0){
		cout<<"ERROR: Partition::Partition(smat,int) : _nitems("<<_nitems<<")!=sitems.size("<<sitems.size()<<")"<<endl;
		exit(1);
	}
}

void Partition::xtrConsPart(multimap<int,string,greaterThan> consPart, int ofs){
	_partitionf=NULL;
	_items_offset=ofs;
	_nclusters=consPart.size();
	_nsingletons=0;
	_nitems=0;
	string str;
	for(multimap<int,string,greaterThan>::iterator c=consPart.begin();c!=consPart.end();c++){
		svect cl;
		stringstream si;
		si<<c->first;
		cl.push_back(si.str());
		stringstream ss;
		ss<<c->second;
		while(ss>>str)cl.push_back(str);
		if(cl.size()-_items_offset==1)_nsingletons++;
		_nitems+=cl.size()-_items_offset;
		clusters.push_back(cl);
	}
		
}

///\f$ H=\sum_k\,p_k\,\log\left(p_k\right) \f$
double Partition::H(){
	double H=0.0;
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++)
		H-=(cl->size()-_items_offset)*log( (double) (cl->size()-_items_offset) );
	return (H/_nitems+log( (double) _nitems ));
}
///\f$ H=\frac{1}{q-1}\,\sum\,p_k^q
double Partition::TS(double q){
//double Partition::TS(double q=EXTENSIVITY_DEFAULT_TSALLIS){
	if(q==1.0)return H();
	double S=0.0;
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++)
		S+=pow( (cl->size()-_items_offset) , q );
	S=1.0-S/pow(_nitems,q)  ;
	return ( S/(q-1.0) ) ;
}
///\f$ H=\frac{1}{q-1}\,\sum\,p_k^q
double Partition::RS(double q){
//double Partition::RS(double q=EXTENSIVITY_DEFAULT_RENYI){
	if(q==1.0)return H();
	double S=0.0;
	///Renyi for q->oo gives -log sup p_k
	if(q>EXTENSIVITY_MAX) { 
		//if(_it_largest_cluster==NULL)cout<<"ERROR: Partition::RS : _it_largest_cluster=NULL"<<endl;
#ifdef DEBUG
		cout<<"DEBUG: Partition::RS : _largest_cluster="<<_largest_cluster<<" _items_offset="<<_items_offset<<" _nitems="<<_nitems<<endl;
		cout<<"DEBUG: max-sup="<<-log( 1.0*(_it_largest_cluster->size()-_items_offset)/_nitems )<<endl;
#endif
		return ( -log( 1.0*(_it_largest_cluster->size()-_items_offset)/_nitems ) );
	}
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++)
		S+=pow( (cl->size()-_items_offset) , q );
	S=log( S/pow(_nitems,q) );
	return ( S/(1.0-q) ) ;
}

double Partition::JQnorm(double q=1.0){
	if(q==1.0)  return exp(H());
	else
		return exp(TS(q));
}

int Partition::getClusterSize(const string& item){
	if(_items_offset!=2){
		cout<<"Error: partitions were defined with no name label : _items_offset="<<_items_offset<<endl;
		exit(1);
	}
	int clsize;
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++){
		clsize=atoi( (*(cl->begin()+0)).c_str() );
		for(svect::iterator it=cl->begin()+_items_offset;it!=cl->end();it++)
			if(strcmp(item.c_str(),it->c_str())==0)return clsize;
	}
	clsize=-1;
	return clsize;
}

string Partition::getClusterName(string& item){
	if(_items_offset!=2){
		cout<<"Error: partitions were defined with no name label : _items_offset="<<_items_offset<<endl;
		exit(1);
	}
	string clname;
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++){
		clname=*(cl->begin()+1);
		for(svect::iterator it=cl->begin()+_items_offset;it!=cl->end();it++)
			if(strcmp(item.c_str(),it->c_str())==0)return clname;
	}
	clname="NAN";
	return clname;
}

int Partition::getClusterIdx(string& item){
	int idx=-1;
	for(int i=0;i<clusters.size();i++)
		for(svect::iterator it=clusters[i].begin();it!=clusters[i].end();it++)
			if(strcmp(item.c_str(),it->c_str())==0)return i;
	return idx;
}

bool Partition::areEquiv(string a, string b){
	if(getClusterIdx(a)==getClusterIdx(b)) return true;
	return false;
}

string Partition::getClusterName(int& clidx){
	if(_items_offset!=2){
		cout<<"Error: partitions were defined with no name label : _items_offset="<<_items_offset<<endl;
		exit(1);
	}
	if(clidx>clusters.size() || clidx<0){
		cout<<"Error : cluster index "<<clidx<<" beyond valid range [0-"<<clusters.size()<<")"<<endl;
		exit(1);
	}
	return clusters[clidx][1];
}

string Partition::areWithinSameCluster(string ita, string itb){
	string clna,clnb;
	clna=getClusterName(ita);
	clnb=getClusterName(itb);
	string noname="NAN";
	if(strcmp(clna.c_str(),clnb.c_str())==0){
		return clna;
	}
	else if(strcmp(clna.c_str(),noname.c_str())==0 ){
		noname="NAN1";
		return noname;
	}
	else if(strcmp(clnb.c_str(),noname.c_str())==0) {
		noname="NAN2";
		return noname;
	}
	return "x";
}

void Partition::edsc(Partition* p2){
	Partition pi=intersection(p2);
	cout<<(2*pi.card()-card()-p2->card() )<<endl;
}

void Partition::td(Partition* p2, double q=1.0){
	cout<<fabs( JQnorm(q)-p2->JQnorm(q) )<<endl;
}
	
void Partition::vipp(Partition* part2){
//double Partition::vipp(Partition& part2){
	int ofs2=part2->cluster_offset();
	long int Nitems=_nitems; //Both partitions are expected to have the same number of items
	long int Nitemsb=part2->_nitems;
	//svect it_found;
	svect it_missing;
	double Iab,tnab;
	double Sa,Sb;
	Sa=Sb=0.0;
	Iab=tnab=0.0;
	_it_found.clear();
	for(smat::iterator cla=clusters.begin();cla!=clusters.end();cla++){				/// For each cluster cla and
		Sa+=( cla->size()-_items_offset )*log( (double) ( cla->size()-_items_offset) );
		for(smat::iterator clb=part2->clusters.begin();clb!=part2->clusters.end();clb++){  	///  clb , find the 
			if(cla==clusters.begin()) Sb+=( clb->size()-ofs2 )*log( (double) ( clb->size()-ofs2) );
			double nab=0;
			for(svect::iterator ita=cla->begin()+_items_offset;ita!=cla->end();ita++){
				for(svect::iterator itb=clb->begin()+ofs2;itb!=clb->end();itb++){
					if(strcmp((*ita).c_str(),(*itb).c_str())==0){
						nab++;							/// number of items they have in common.
						_it_found.push_back(*ita);
					}
					//nab+=strcmp((*ita).c_str(),(*itb).c_str())==0?1.0:0.0;
				}
			}
			tnab+=nab;
			if(nab>0) Iab+=nab*log(nab);
		}
	}
	//missing(&_it_found);
	Sa=( Sa/Nitems - log( (double) Nitems ) );
	if(Nitems==Nitemsb){					/// VI make sense only if number items is equeal.
		//Sb=( Sb/Nitemsb - log( (double) Nitemsb ) ); 	
		Sb=( Sb - Nitemsb*log( (double) Nitems ) )/Nitems; 	/// Using the same normalization for both entropies
	} else {
		Sb=( Sb - Nitemsb*log( (double) Nitems ) )/Nitems; 	/// Using the same normalization for both entropies
		tnab+=Nitems-Nitemsb;                          /// is tantamout to assuming each item from partition1 not present in partition2 is a singleton in partition2
		cout<<"#Nitems= "<<Nitems<<" <> Nitemsb="<<Nitemsb<<" tnab="<<tnab<<endl;
	}
	Iab=( Iab-tnab*log((double)Nitems) )/(double)Nitems;	// Iab = Sum_cla Sum_clb nab/N * log( nab/N )
	cout<<(Sa + Sb - 2.0*Iab )<<"\t"<<Sa<<"\t"<<Sb<<"\t"<<Iab<<endl;
//	return (Sa + Sb - 2.0*Iab ) ;
}

void Partition::missing(){
	missing(&_it_found);
}
void Partition::missing(svect* it_found){
	if(it_found->size()==0) return;
	svect it_missing;
        svect::iterator msit;
        for(smat::iterator cla=clusters.begin();cla!=clusters.end();cla++){
                for(svect::iterator ita=cla->begin()+_items_offset;ita!=cla->end();ita++){
                        it_missing.push_back(*ita);
                        msit=it_missing.end();
                        for(svect::iterator itb=it_found->begin();itb!=it_found->end();itb++)
                                if(strcmp((*ita).c_str(),(*itb).c_str())==0)
                                        it_missing.erase(msit);
                }
	}
	if(it_missing.size()>0){
		cout<<"#Partition1 items MISSING in partition2 : ";
		for(svect::iterator it=it_missing.begin();it!=it_missing.end();it++)
			cout<<*it<<" ";
		cout<<endl;
	}
}

void Partition::SubsProject(sset& itemset){
	if(VERBOSE)cout<<"#Erase : ";
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++)
		for(svect::iterator it=cl->begin()+_items_offset;it!=cl->end();it++)
			if(itemset.find(*it)==itemset.end()){ //If item it isn't in itemset
				if(VERBOSE)cout<<*it<<" ";
				it=cl->erase(it);		//erase it.
				_nitems--;
				--it;				//Then adjust appropriately the iterator to point to previous item in cl,
				int i=atoi( (*cl)[0].c_str() ); //get the number of elements in cl (assumes partFmt=PART) and
				stringstream ss;
				ss<<--i;
				if(i>0)(*cl)[0]=ss.str();		//decreas it by one if there are more elements remaining, otherwise
				else {
					cl=clusters.erase(cl)-1; //delete cluster (it was a singleton) and jump to next cluster.
					if(VERBOSE)cout<<"(Singleton) ";
					break;
				}
				if(VERBOSE)cout<<"(New cluster size= "<<(*cl)[0]<<") ";
			}
	if(VERBOSE)cout<<endl;
	_resetMembers();
}

Partition Partition::intersection(Partition* part2){
	int ofs2=part2->cluster_offset();
	int cln=0;
	long int nitems=part2->n_items();
	if(nitems<_nitems)nitems=_nitems;
	string isize,iname; 
	smat intersection;
	for(smat::iterator cla=clusters.begin();cla!=clusters.end();cla++){
		for(smat::iterator clb=part2->clusters.begin();clb!=part2->clusters.end();clb++){
			svect inter_cl; /// Will contain intersection between cla and clb
			for(svect::iterator ita=cla->begin()+_items_offset;ita!=cla->end();ita++)
				for(svect::iterator itb=clb->begin()+ofs2;itb!=clb->end();itb++)
					if(strcmp((*ita).c_str(),(*itb).c_str())==0)
						inter_cl.push_back(*ita);
			if(inter_cl.size()>0){ /// Set its size and name. These will be its first two entries
				cln++;
				stringstream ssn;
				ssn<<"N"<<cln;	/// First, the name of the cluster
				iname=ssn.str();	/// the name of the cluster
				svect::iterator it;
				if(_items_offset>1||ofs2>1){
					it=inter_cl.begin();
					inter_cl.insert(it,1,iname); 
				}
				stringstream sss;/// Second, the size of the cluster
				int i=(int)(inter_cl.size()-_items_offset+1);
				sss<<i;
				isize=sss.str();/// Second, the size of the cluster
				it=inter_cl.begin();
				inter_cl.insert(it,1,isize);
				stringstream rsss; ///Third, the fake reverse-sorting label
				rsss<<nitems-i;
				isize=rsss.str();
				it=inter_cl.begin();
				inter_cl.insert(it,1,isize); ///Third, the fake reverse-sorting label
			}
			if(inter_cl.size()>0)intersection.push_back(inter_cl);
		}
	}
	if(intersection.size()>0){	///Reverse sort : Larger clusters firstl singletons last
		sort(intersection.begin(),intersection.end());
		int maxitems=0;
		for(smat::iterator cl=intersection.begin();cl!=intersection.end();cl++){
			cl->erase(cl->begin()); ///REmove the first fake reverse-sorting label of each cluster
			if(cl->size()-_items_offset>maxitems){
				//maxitems=sitems.size();
				maxitems=cl->size()-_items_offset;
				_largest_cluster=clusters.size()-1; //index refering to largest cluster
				_it_largest_cluster=clusters.begin(); ///pointer refering to largest cluster
			}
		}
	}
	return Partition(&intersection, ofs2>_items_offset?ofs2:_items_offset);
}

Partition Partition::operator+(Partition& part2){
	int ofs2=part2.cluster_offset();
	int cln=0;
	long int nitems=part2.n_items();
	if(nitems<_nitems)nitems=_nitems;
	string isize,iname; 
	smat enose;
	map<long int, smat > pan_enose;
	svect::iterator iab,iae,ibb,ibe;
	iab=clusters[largest_Cluster()].begin()+_items_offset; ///Start of largest cluster in this
	iae=clusters[largest_Cluster()].end();			///and end.
	svect acl (iab,iae);	///extract only the elements, i.e., w/o the labels.
	ibb=part2.clusters[part2.largest_Cluster()].begin()+ofs2; ///Start of largest cluster in [art2
	ibe=part2.clusters[part2.largest_Cluster()].end();	///and its end.
	svect bcl (ibb,ibe);	///extract only the elements, i.e., w/o the labels.
	svect enose_cl=acl+bcl; ///First join the two largest clusters
	enose.push_back(enose_cl); 
	if(enose_cl.size()<nitems)enose.push_back(enose_cl);
	for(smat::iterator cla=clusters.begin();cla!=clusters.end();cla++){
		if((*cla)<=enose_cl)continue;
		for(smat::iterator clb=part2.clusters.begin();clb!=part2.clusters.end();clb++){
			if((*clb)<=enose_cl)continue;
			svect sum_cl=(*cla)+(*clb);
			enose.push_back(sum_cl); 
		}	
	}
}

bool isInteger(string& str){
	for(int i=0;i<str.length();i++)
		if(!isdigit(str.at(i)))return false;
	return true;
}

void Partition::mclTabFile(char* mcltabf){
	_mcltabf=mcltabf;
	if(!QUIET)cout<<"#Reading MCL tab file "<<string(_mcltabf)<<endl;
	ifstream is(_mcltabf);
	if(!is){
		cout<<"ERROR: Partition::mclTabFile : Cannot open file "<<string(_mcltabf)<<endl;	
		exit(1);
	}
	string it;
	int c=0;
	while(is>>it){
		//if(!isdigit(it.at(0))){
		if(!isInteger(it)){
			getline(is,it);
			continue;
		}
		else{
			int idx=atoi(it.c_str());
			is>>it;
			_mcltab[idx]=it;
			c++;
		}
	}
	if(n_items()>0 && c!=n_items()){
		cout<<"ERROR: Partition::mclTabFile : number of keys-items pairs doesn't match number of items read for partition."<<endl;
		exit(1);
	}
	if(!QUIET)cout<<"#MCL tab file: found "<<c<<" key-label pairs. Last label is "<<_mcltab[c-1]<<endl;
}
void Partition::swapLabels(char* mcltabf=NULL){
	if(mcltabf)mclTabFile(mcltabf);
	if(_mcltabf==NULL){
		cout<<"#WARNING: Couldn't translate items labels. NO tab file defined."<<endl;
		return ;
	}
	smat nclusters;
	smat* cls=&clusters;
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++){
		svect nits;
		for(svect::iterator it=cl->begin();it!=cl->end();it++){
			if(it>=cl->begin()+_items_offset) 
				nits.push_back( _mcltab[ atoi( (*it).c_str() ) ]  );
			else nits.push_back(*it);
		}
		nclusters.push_back(nits);
	}
	//if(!QUIET)cout<<"#Items translated according to tabfile. Previous Last item= "<<*((clusters.rbegin())->rbegin())<<endl;
	clusters=nclusters;
	_resetMembers();
	if(!QUIET)cout<<"#Items translated according to tabfile. New Last item= "<<*((nclusters.rbegin())->rbegin())<<endl;
}

void Partition::partitionInputFormat(partFileFormat iformat=partFmtPART){
	switch(iformat){
		case partFmtMCL:
			if(VERBOSE)cout<<"#Reading partition input format MCL"<<endl;
			break;
		case partFmtFREE:
			if(VERBOSE)cout<<"#Reading partition input format FREE"<<endl;
			break;
		case partFmtPART:
		default:
			if(VERBOSE)cout<<"#Reading partition input format PARTANALYZER"<<endl;
			break;
	}
}

const char* reol_="$\0";
//bool _readStringf(ifstream& _is, string& it, partFileFormat& iformat, int nit=1){ //Default input format is partanalyzer's own format
bool _readStringf(ifstream& _is, string& it, partFileFormat& iformat, int nit=1, const char* eol=reol_){ //Default input format is partanalyzer's own format
	if(!(_is>>it)) return false;
	if(it.compare("(mclheader")==0){
		iformat=partFmtMCL;
		if(!QUIET)cout<<"#Detected partition input format MCL"<<endl;
	}
	switch(iformat){
		case partFmtMCL:
			if( !isdigit(it.at(0)) ){
#ifdef DEBUG
				cout<<"#DEBUG: not a digit seen {"<<it<<"}: nit="<<nit<<endl;
#endif
				string line;
				getline(_is,line);
				if(_readStringf(_is,it,iformat,nit) ) 
					return true;
				else return false;
			}
			//if( isdigit(it.at(0)) ){
			else if(nit>0) return true;
			else{
#ifdef DEBUG
				cout<<"#DEBUG: digit seen {"<<it<<"}: nit="<<nit<<endl;
#endif
				string line;
				int cispos=_is.tellg();
				cispos-=it.length();
				//getline(_is,line,'$');
				getline(_is,line,eol[0]);
				_is.seekg(cispos);
				stringstream ss;
				ss<<line;
#ifdef DEBUG
				cout<<"#DEBUG: line seen : {"<<line<<"}"<<endl;
				cout<<"#DEBUG: stream seen : {"<<ss.str()<<"}"<<endl;
#endif
				line.clear();
				long int nitems=0;
#ifdef DEBUG			
				cout<<"#DEBUG: counting clusters elements... ";
#endif
				//while(ss>>line && line.compare("$")!=0) ///In MCL, the first number reefers to clusters name; line ends with $
				while(ss>>line && line.compare(eol)!=0) ///In MCL, the first number reefers to clusters name; line ends with $
				{
#ifdef DEBUG
					cout<<line<<" ";
#endif
					nitems++;	
				}
#ifdef DEBUG
				cout<<"\n#DEBUG: nitems count : "<<nitems<<endl;
#endif
				stringstream iss;
				iss<<nitems;
				it=iss.str();
				return true;
			}
			break;
		case partFmtFREE:
			cout<<"ERROR: input format FREE not yet implemented"<<endl;
			exit(1);
			break;
			/*
			*/
		case partFmtPART:
		default:
			if(strcmp(it.substr(0,1).c_str(),"#")==0){
				if(VERBOSE)cout<<"#Skipping comment... \n"<<it<<" ";
				string line;
				getline(_is,line);
				//ostringstream oss;
				//while(oss<<it)
				//	if(VERBOSE)cout<<oss.str()<<" ";
				if(VERBOSE)cout<<endl;
				if(_readStringf(_is,it,iformat) ) 
					return true;
				else return false;
			}
			else if(iformat==partFmtFREE){
#ifdef DEBUG
				cout<<"#DEBUG: element seen {"<<it<<"}: nit="<<nit<<endl;
#endif
				string line;
				int cispos=_is.tellg();
				cispos-=it.length();
				getline(_is,line);
				_is.seekg(cispos);
#ifdef DEBUG
				cout<<"#DEBUG: line seen : {"<<line<<"}"<<endl;
#endif
#ifdef DEBUG			
				cout<<"#DEBUG: counting clusters elements... ";
#endif
				long int nitems=1;
				int i=-1;
				while(line.at(++i)) ///In FREE, simply elements are listed.
				{
#ifdef DEBUG
					cout<<line.at(i);
#endif
					if(isspace(line.at(i))) nitems++;	
				}
#ifdef DEBUG
				cout<<"\n#DEBUG: nitems count : "<<nitems<<endl;
#endif
				stringstream iss;
				iss<<nitems;
				it=iss.str();
				return true;
				
			}
			break;
	}
	return true;
}
bool _isComment(ifstream& _is, string& it, partFileFormat& iformat, int nit, char* eol="$\0"){
	switch(iformat){
		case partFmtMCL:
			if( !isdigit(it.at(0) ) ){
				string line;
				getline(_is,line);
				if(_readStringf(_is,it,iformat,nit) ) 
					return true;
				else return false;
			}
			break;
		case partFmtPART:
			break;
	}
}

void Partition::_readClusters(){ ///For the time being, we'll assume each cluster has its number of items as the first string, its name as the second and then the items we'll follow:
	string it;		///That is, e.g, "95 N245 Alx3 Alx4 Cart1 ..." This would require an offset _items_offset=2. However, the lowest offset we expect is 
	int nit,nitems;	/// _items_offset=1, as at least the first item must be the cluster size. Including the name is optional.
	_nsingletons=0;
	_npairs=0;
	int maxitems=0;
	svect _items;
	ifstream _is(_partitionf);
        if(!_is){
                cout<<"Cannot open file "<<string(_partitionf)<<endl;
                exit(1);
        }
	nit=0; ///Item-read number
	_nitems=0;
	if(!QUIET) cout<<"#Reading clusters from partition "<<_partitionf<<endl;
	partitionInputFormat(_piformat);

	//while(_is>>it && strcmp(it.substr(0,1).c_str(),"#")!=0){ //Comments only at the end of file and starting with #
	while(_readStringf(_is,it,_piformat,nit) ){ //Comment anywhere starting with #
#ifdef DEBUG
			cout<<"Seen: ("<<it.substr(0,1)<<")";
#endif
			if(nit==0){
				nitems=atoi(it.c_str());	///Number of items to be read.
				_items.assign(nitems + _items_offset,""); ///create items list accordingly, leaving space for the additional (_items_offset) indexes.
				_items[nit]=it;				///save the first item
#ifdef DEBUG
				cout<<" (_items.size="<<_items.size()<<")"<<"n="<<it<<"\t";
#endif
			}
			//while(nit<nitems+_items_offset-1 && _is>>it){ ///If we haven't read all records, keep reading and setting up the cluster
			while(nit<nitems+_items_offset-1 && _readStringf(_is,it,_piformat) ){ ///If we haven't read all records, keep reading and setting up the cluster
				nit++;
				_items[nit]=it;
				if(nit>_items_offset-1){
					sitems.insert(it);///Update set ot items with new element.
					_vitems.push_back(it);
#ifdef DEBUG
					cout<<"[nit="<<nit<<"/sitems<<"<<it<<"]";
#endif
				}
#ifdef DEBUG
				cout<<it<<"\t";
#endif
			}
			if(nit<nitems+_items_offset-1){
				cout<<"#"<<program<<": Reading cluster: last item read: "<<_items[nit]<<". Skipping this line and rest of input file."<<endl;
				break;
			}
#ifdef DEBUG
			cout<<" length="<<_items.size();
			cout<<endl;
#endif
			if(_items.size()==1+_items_offset){
				_nsingletons++;
				ssingletons.insert(_items[_items_offset]);
			}
			if(_items.size()==2+_items_offset)_npairs++;
			clusters.push_back(_items);
			_nitems+=_items.size()-_items_offset;
			nit=0;
			if(_items.size()>maxitems){
				maxitems=_items.size();
				_largest_cluster=clusters.size()-1; //index refering to largest cluster
				_it_largest_cluster=clusters.begin(); ///pointer refering to largest cluster
			}
	}
	if(VERBOSE) cout<<"#Finish Reading partition: "<<endl;
	if(!QUIET) cout<<"#items= "<<_nitems<<" clusters= "<<clusters.size()<<" #singletons= "<<_nsingletons<<" #pairs= "<<_npairs<<" \%non-trivial= "<<(_nsingletons+_npairs)*1.0/clusters.size()<<" largest-cluster-size= "<<maxitems-_items_offset<<" largest-cluster-index= "<<_largest_cluster<<" last-cluster-size= "<<_items.size()-_items_offset<<" Last-item-read= "<<_items[_items.size()-1]<<endl;
	_nclusters=clusters.size();
	if(_nsingletons!=ssingletons.size()){
		cout<<"ERROR: _readClusters : Not a sound Partition : non mutually disjoint singleton clusters"<<endl;
		cout<<ssingletons<<endl;
		exit(1);
	}
}

void Partition::printPartition(){
	partFileFormat format=_piformat;
	printPartition(format);
}

void Partition::printPartition(partFileFormat format){
	int ofs=0;
	switch(format){
		case partFmtMCL:
			{
				if(!QUIET)cout<<"#BeginMCLPartitionMatrix"<<endl;
				ostringstream oss;
				cout<<"(mclheader\nmcltype matrix\ndimensions "<<n_items()<<"x"<<n_clusters()<<"\n)"<<endl;
				cout<<"(mclmatrix\nbegin"<<endl;
				long int it=0;
				for(int i=0;i<clusters.size();i++){
					cout<<i<<"\t";
					for(int j=_items_offset;j<clusters[i].size();j++,it++)
						if(_mcltabf!=NULL ){
							for(map<int,string>::iterator iter=_mcltab.begin();iter!=_mcltab.end();iter++)
								if((*iter).second.compare(clusters[i][j])==0)cout<<(*iter).first<<" ";
						}else{
							cout<<it<<" ";
							oss<<it<<"\t"<<clusters[i][j]<<endl; //If mcltab file was not declared, built one ad-hoc
						}
					cout<<"$"<<endl;
				}
				cout<<")"<<endl;
				if(!QUIET)cout<<"#EndMCLPartitionMatrix"<<endl;
				//if(_mcltabf!=NULL ){
					if(!QUIET)cout<<"#BeginMCLTabfile"<<endl;
					cout<<"#Tab file of partition "<<(FileName()==NULL?"":FileName())<<endl;
					cout<<oss.str();
					oss.flush();
					if(!QUIET)cout<<"#EndMCLTabfile"<<endl;
				//}
				break;
			}
		case partFmtFREE:
			ofs=_items_offset;
		case partFmtPART:
		default:
			if(!QUIET)cout<<"#Clusters: "<<clusters.size()<<endl;
			for(int i=0;i<clusters.size();i++){
				for(int j=0+ofs;j<clusters[i].size();j++)
					//cout<<"("<<i<<","<<j<<")="<<clusters[i][j]<<"\t";
					if(_mcltabf!=NULL && j>=_items_offset)
						cout<<_mcltab[ atoi( clusters[i][j].c_str() ) ]<<"\t";
					else
						cout<<clusters[i][j]<<"\t";
				cout<<endl;
			}
			break;
	}
}

//neighborhood Partition::Neighborhood(string item)
svect Partition::getClusterOf(string item)
{
	//sset* mycluster=NULL;
	//double nf=1.0;
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++)
		for(svect::iterator it=cl->begin()+_items_offset;it!=cl->end();it++)
			if(item.compare(*it)==0) {
					svect mycluster (cl->begin()+_items_offset,cl->end());
					return mycluster;
			}
	/*
			{
				for(it=cl->begin();it!=cl->end();it++)
					mycluster->insert(*it);
				break;
			}
	if(!mycluster){
		mycluster->insert("NAN");
		nf=0.0;
	}
	return pair<string,pair<sset,double> > (item,pair<sset,double>(*mycluster,nf));
	*/
	svect mycluster;
	return mycluster;
}

sset Partition::getItems(cluster* cl)
{
	sset clItems;
	for(cluster::iterator it=cl->begin();it!=cl->end();it++)
		clItems.insert(*it);
	return clItems;
}

sset Partition::getItems()
{
	sitems.clear();
	for(smat::iterator cl=clusters.begin();cl!=clusters.end();cl++)
		for(svect::iterator it=cl->begin()+_items_offset;it!=cl->end();it++)
			sitems.insert(*it);
#ifdef DEBUG
	cout<<"Partition::getItems : Items partition ";
	if(FileName()!=NULL)cout<<FileName();
	cout<<" : "<<endl;
	cout<<"{";
	for(sset::iterator it=sitems.begin();it!=sitems.end();it++)
		cout<<*it<<",";
	cout<<"}"<<endl;
#endif
	return sitems;
}

bool Partition::operator==(Partition& part2)
{
	bool found;
	if(_nclusters!=part2.n_clusters())return false;
	if(_nsingletons!=part2.n_singletons()) return false;
	if(_nitems!=part2.n_items()) return false;
	for(smat::iterator cla=clusters.begin();cla!=clusters.end();cla++)
		for(svect::iterator ita=cla->begin()+_items_offset;ita!=cla->end();ita++){
			found=false;
			for(smat::iterator clb=part2.clusters.begin();clb!=part2.clusters.end();clb++)
				for(svect::iterator itb=clb->begin()+_items_offset;itb!=clb->end();itb++)
					if(ita->compare(*itb)==0) found=true;
			if(!found)return false;
		}
	return true;
}

bool Partition::operator<=(Partition& part2){
	if(n_clusters()<part2.n_clusters())return false;
	//if(*this==part2)return true;
	long int clfound=0;
	for(smat::iterator cla=clusters.begin();cla!=clusters.end();cla++)
		for(smat::iterator clb=part2.clusters.begin();clb!=part2.clusters.end();clb++){
			svect acl (cla->begin()+_items_offset,cla->end());
			svect bcl (clb->begin()+part2.cluster_offset(),clb->end());
			//if(*cla<=*clb)clfound++;
			if(acl<=bcl)clfound++;
		}
	if(clfound==n_clusters())return true;
	return false;
}

/*Partition Partition::operator*(Partition& part2){
	return intersection(&part2);
}
*/
vector<double > Partition::purityScore(Partition* part2){
	double strict=0.0;
	double lax=0.0;
	long int itfound;
	vector<double > pscores;
	for(smat::iterator cla=clusters.begin();cla!=clusters.end();cla++){
		svect acl (cla->begin()+_items_offset,cla->end());
		if(acl.size()==1)continue;
		for(smat::iterator clb=part2->clusters.begin();clb!=part2->clusters.end();clb++){
			svect bcl (clb->begin()+part2->cluster_offset(),clb->end());
			if(bcl.size()==1)continue;
			if(acl<=bcl)lax++;
			if(acl==bcl)strict++;
		}
	}
	pscores.push_back(strict/(double)part2->n_nonSingClusters()); ///Normalize by the reference partition
	pscores.push_back(lax/(double)n_nonSingClusters()); ///Normalize by the target partition
	return pscores;
}

bool Partition::isaPartition(){
	svect its;
	if(sitems.size()==0)getItems();
	for(sset::iterator it=sitems.begin();it!=sitems.end();it++)
		its.push_back(*it);
	return isaPartitionOf(its);
}
bool Partition::isaPartitionOf(svect& svecOfelements){///One single clusters does make a sound partitions
	svect unionOfcls;
	if(clusters.size()==1){
		//return true;  ///We shouldn't need to go beyond this, if we can make sure sitems comes always well defined.
		unionOfcls=clusters[0];	
		for(int i=0;i<_items_offset;i++)unionOfcls.erase(unionOfcls.begin());
	}else for(int i=0;i<clusters.size()-1;i++){
	///With the introduction of Partition::ssingletons, no instantiated Partition can have the same singleton element more than once.
	///Hence, there is no need in checking singletons againts singletons.
	///Looping over them means that, on an Opteron 846 2GH, it takes 15sec to check a partition with 502 elements and
	///on average 185 clusters and 95 singletons. Does skipping singletons helps to significantly improve this figure?
	///Also, take svect cla= clusters[i]; out of the second loop (that was weird!).
	///Ok, preliminary test shows an increase in speed of ~28%, i.e, 9.5sec/partition now versus 13.2sec before! (Compared both ways using 
	///a set of 19 partitions of those same sizes.
		svect cla= clusters[i];
		///For the moment, lets assume clusters are not ordered in descending order by size
		if(cla.size()==1+_items_offset)continue;///SKIP check for singletons, but non-singletons ARE checked againts clb singletons!!
		//if(cla.size()==1+_items_offset)break;///Makes no difference with respect to 'continue' for a set of 19 partitions of those sizes
		for(int j=i+1;j<clusters.size();j++){
			svect clb= clusters[j];
			int n=0;
			while(_items_offset>n++) {
				if(j==i+1)cla.erase(cla.begin()); //Use the same loop as for clb, but got to do it only once for cla!!
				clb.erase(clb.begin()); 
			}
//	for(smat::iterator cla=clusters.begin();cla!=clusters.end();cla++)
//		for(smat::iterator clb=clusters.begin();clb!=clusters.end();clb++){
//			if(cla==clb)continue;
			svect icl=cla*clb;
			if(icl.size()>0){
				cout<<"#ERROR: isaPartitionOf : Non mutually disjoint clusters found : "<<i<<"("<<cla[0]<<" ...) - "<<j<<"("<<clb[0]<<" ...)"<<endl;
				cout<<"cla="<<cla<<endl;
				cout<<"clb="<<clb<<endl;
				cout<<"cli="<<icl<<endl;
				return false;
			}
			svect sum=(cla+clb); 
			unionOfcls+=sum;
		}
	}
	if(n_singletons()!=n_clusters() && !(unionOfcls==svecOfelements)){
		cout<<"#ERROR: isaPartitionOf : Clusters' union ("<<unionOfcls.size()<<") does  not cover the whole underlying set of elements ("<<svecOfelements.size()<<")"<<endl;
		cout<<"#Union="<<unionOfcls<<"\n#Set="<<svecOfelements<<endl;	
		return false;
	}
	return true;
}

#endif //END  _CLASS_PARTITION

