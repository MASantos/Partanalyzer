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


#ifndef _CLASS_PARTITIONSTATS
#define _CLASS_PARTITIONSTATS 1

#include "PartitionStats.h"

PartitionStats::PartitionStats(vector<Charr > fnames, partFileFormat iformat=partFmtPART, double extensivity=EXTENSIVITY_DEFAULT, int ofs=2,int clstat_normalization_ofs=0){
	_clsnofs=clstat_normalization_ofs;
	_fnamel=fnames;
	_npart=0;
	_DIST_SUBSPROJECT=false;
	_entensivity_degree=extensivity;
	if(VERBOSE){
		cout<<"#Instantiating each partition ("<<fnames.size()<<"/"<<_fnamel.size()<<") ..."<<endl;
		cout<<"#First partition... "<<fnames[0].car<<"/"<<_fnamel[0].car<<endl;
		cout<<"#Last partition... "<<fnames[fnames.size()-1].car<<"/"<<_fnamel[_fnamel.size()-1].car<<endl;
	}
	
	for(vector<Charr >::iterator f=_fnamel.begin();f!=_fnamel.end();f++){
		_npart++;
		if(!QUIET)cout<<"#Partition "<<_npart<<endl;
		Partition q=Partition(f->car,iformat,ofs);
		_partitionl.push_back(q);
		_hasseNodes[q.n_clusters()].push_back(q);
	}
	if(_npart<1){
		cout<<"ERROR: invalid number of partitions "<<_npart<<"<1 : Check value provided for option norm"<<endl;
		exit(1);
	}
	if(_npart>300&&!QUIET){
		cout<<"#Finished instantiating all partitions."<<endl;
		systemDate();
	}
}

///Distance metric induced by the cardinality : \f$d\left(P,P'\right)\,\equiv\,2\,card\left(P,P'\right)\,-\,card\left(P\right)\,-\,card\left(P'\right)\f$
long int PartitionStats::_pmetric(Partition& p1, Partition& p2, const int f){
	return ( 2*card(p1*p2)-card(p1)-card(p2) );
}
/* ///Distance metric induced by the entropy : \f$d\left(P,P'\right)\,\equiv\,2\,H\left(P,P'\right)\,-\,H\left(P\right)\,-\,H\left(P'\right)\f$
double PartitionStats::_pmetric(Partition& p1, Partition& p2, const double f=0.0){
	return ( 2.0*H(p1*p2)-H(p1)-H(p2) );
}*/
///Distance metric induced by the entropy : \f$d\left(P,P'\right)\,\equiv\,2\,H\left(P,P'\right)\,-\,H\left(P\right)\,-\,H\left(P'\right)\f$
double PartitionStats::_pmetric(Partition& p1, Partition& p2, pmetricv metric){
	if(metric==jeffreyQnorm) return fabs(log(_f(p1,metric))-log(_f(p2,metric))); //Tarantola distance
	Partition q=p1*p2;
	return ( 2.0*_f(q,metric)-_f(p1,metric)-_f(p2,metric) );
	//return ( 2.0*_f(p1*p2,metric)-_f(p1,metric)-_f(p2,metric) );
}
double PartitionStats::_f(Partition& p, pmetricv metric){ //POTENTIAL ASSOCIATED WITH A PARTITION.
	double f;
	switch(metric){
		case boltzmann:
				f=BK(p);
				break;
		//case entropy:
		case shannon:
				f=H(p);
				break;
		case cardinality:
				if(VERBOSE)cout<<"#WARNING: PartitionStats::_f : using cardinality as double"<<endl;
				f=card(p);
				break;
		case tsallis:	
				if(_entensivity_degree==EXTENSIVITY_DEFAULT)
					f=TS(p);
				else
					f=TS(p,_entensivity_degree);
				break;
		case renyi:
				if(_entensivity_degree==EXTENSIVITY_DEFAULT)
					f=RS(p);
				else
					f=RS(p,_entensivity_degree);
				break;
		case jeffreyQnorm:
				f=JQnorm(p,_entensivity_degree);
				break;
		default:
				cout<<"ERROR: PartitionStats::_f : unknown metric."<<endl;
				exit(1);
				break;
	}			
	return f;
}

void PartitionStats::iPotential(pmetricv metric=renyi){
	if(!QUIET)cout<<"#BeginIPotential"<<endl;
	for(vector<Partition >::iterator p=_partitionl.begin();p!=_partitionl.end();p++)
		cout<<p->FileName()<<"\t"<<_f(*p,metric)<<endl;
	if(!QUIET)cout<<"#EndIPotential"<<endl;
}

void PartitionStats::pmeasures(pmeasure measure , pmetricv metric){
	distancesPrintHeadComment(measure, metric,(flagheader) BEGIN);
	for(int i=0;i<_partitionl.size()-1;i++){
		Partition pa=_partitionl[i];
		for(int j=i+1;j<_partitionl.size();j++){
			if(_DIST_SUBSPROJECT) pa=_partitionl[i];
			Partition pb=_partitionl[j];
			Partition pi=pa*pb;
			if(_DIST_SUBSPROJECT){
				if(VERBOSE)cout<<"#Projecting onto common subspace ("<<pi.sitems.size()<<") "<<endl;
				pa.SubsProject(pi.sitems);
				pb.SubsProject(pi.sitems);
			}
			double jointH=_f(pi,metric);
			double ipota=_f(pa,metric);
			double ipotb=_f(pb,metric);
			double vmeasure_h=ipotb<1.0E-12?1.0:1-(jointH-ipota)/ipotb;
			double vmeasure_c=ipota<1.0E-12?1.0:1-(jointH-ipotb)/ipota;
			switch(measure){
				case conditionalEntropy:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<jointH-_f(pb,metric)<<"\t"<<jointH-_f(pa,metric)<<"\t"<<2*jointH-_f(pb,metric)-_f(pa,metric)<<endl;
					break;
				case jointEntropy:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<jointH<<endl;
					break;
				case vmeasureArithmetic:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<0.5*(vmeasure_h+vmeasure_c)<<endl;
					break;
				case vmeasureGeometric:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<sqrt(vmeasure_h*vmeasure_c)<<endl;
					break;
				case vmeasureHarmonic:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<(1.0+beta)*vmeasure_h*vmeasure_c/(beta*vmeasure_h+vmeasure_c)<<endl;
					break;
				case symmetricPurity:{
					double pl=0.0;
					double ps=0.0;
					vector<double> pscores;
					pscores=pa.purityScore(&pb);
					ps+=pscores[0];
					pl+=pscores[1];
					pscores=pb.purityScore(&pa);
					ps+=pscores[0];
					pl+=pscores[1];
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<0.5*ps<<"\t"<<0.5*pl<<endl;
					}
					break;
				default:
					cout<<"ERROR: PartitionStat::measure() : unknown measure case"<<endl;
					exit(1);
			}
		}
	}	
	distancesPrintHeadComment(measure, metric,(flagheader) END);
	_DIST_SUBSPROJECT=false;
}

void PartitionStats::pmeasuresRef(pmeasure measure , pmetricv metric){
	Partition p=_partitionl[0];
	bool usingREF=true;
	distancesPrintHeadComment(measure, metric,(flagheader) BEGIN);
	for(int i=1;i<_partitionl.size()-1;i++){
		for(int j=i+1;j<_partitionl.size();j++){
			if(_DIST_SUBSPROJECT) p=_partitionl[i];
			Partition pb=_partitionl[j];
			Partition pi=p*pb;
			if(_DIST_SUBSPROJECT){
				if(VERBOSE)cout<<"#Projecting onto common subspace ("<<pi.sitems.size()<<") "<<endl;
				p.SubsProject(pi.sitems);
				pb.SubsProject(pi.sitems);
			}
			double jointH=_f(pi,metric);
			double ipota=_f(p,metric);
			double ipotb=_f(pb,metric);
			double vmeasure_h=ipotb<1.0E-12?1.0:1-(jointH-ipota)/ipotb;
			double vmeasure_c=ipota<1.0E-12?1.0:1-(jointH-ipotb)/ipota;
			switch(measure){
				case conditionalEntropy:
					cout<<pb.FileName()<<"\t"<<jointH-_f(p,metric)<<"\t"<<2.0*jointH-_f(p,metric)-_f(pb,metric)<<endl;
					break;
				case jointEntropy:
					cout<<pb.FileName()<<"\t"<<jointH<<endl;
					break;
				case vmeasureArithmetic:
					cout<<pb.FileName()<<"\t"<<0.5*(vmeasure_h+vmeasure_c)<<endl;
					break;
				case vmeasureGeometric:
					cout<<pb.FileName()<<"\t"<<sqrt(vmeasure_h*vmeasure_c)<<endl;
					break;
				case vmeasureHarmonic:
					cout<<pb.FileName()<<"\t"<<(1.0+beta)*vmeasure_h*vmeasure_c/(beta*vmeasure_h+vmeasure_c)<<endl;
					break;
				case symmetricPurity:{
					double pl=0.0;
					double ps=0.0;
					vector<double> pscores;
					pscores=p.purityScore(&pb);
					ps+=pscores[0];
					pl+=pscores[1];
					pscores=pb.purityScore(&p);
					ps+=pscores[0];
					pl+=pscores[1];
					cout<<pb.FileName()<<"\t"<<0.5*ps<<"\t"<<0.5*pl<<endl;
					}
					break;
				default:
					cout<<"ERROR: PartitionStat::measureRef() : unknown measure case"<<endl;
					exit(1);
			}
		}
	}	
	distancesPrintHeadComment(measure, metric,(flagheader) END);
	_DIST_SUBSPROJECT=false;
}

int PartitionStats::arePartitions(){
	int nopart=0;
	int np=1;
	for(vector<Partition >::iterator p=_partitionl.begin();p!=_partitionl.end();p++,np++)
		if(!p->isaPartition()){
			nopart++;
			cout<<"ERROR:\t"<<np<<"\t"<<p->FileName()<<endl;
		}else if(!QUIET) cout<<"#Partition "<<np<<" OK"<<endl;
	return nopart;
}

void PartitionStats::getCover(){
	if(VERBOSE)
		cout<<"#Counting neighborhoods..."<<endl;
	for(vector<Partition >::iterator p=_partitionl.begin();p!=_partitionl.end();p++){
		sset its=p->getItems();
		for(sset::iterator it=its.begin();it!=its.end();it++){
			svect cl=p->getClusterOf(*it);
#ifdef DEBUG
			cout<<"#DEBUG: cluster for item "<<*it<<" : ";
			for(svect::iterator i=cl.begin();i!=cl.end();i++)
				cout<<*i<<" ";
			cout<<endl;
#endif
			sNeighborhood itemNeighborhood(*it,&cl);
			_cover[*it]+=itemNeighborhood;
			if(VERBOSE){
				cout<<"#Current map of Neighborhoods "<<*it<<" : ";
				_cover[*it].gNeighborsList();
			}
		}
	}
}

void PartitionStats::printCover(bool PRINTCONSENSUP=false){
	long int npart=_npart+_clsnofs;
	cout<<"#BeginCover"<<endl;
	cout<<"#Cover size= "<<_cover.size()<<endl;
	map<string,int> ccl;
	int cclso=0;
	for(map<string,sNeighborhood>::iterator it=_cover.begin();it!=_cover.end();it++){
		cout<<it->first;
		cout<<"\t"<<(it->second).eps<<"\t"<<(it->second).eps/(double)npart*100.0;
		cout<<"\t"<<(it->second).neighbors.size();
		stringstream ss;
		for(sset::iterator p=(it->second).neighbors.begin();p!=(it->second).neighbors.end();p++){
			cout<<"\t"<<*p;
			if(PRINTCONSENSUP)ss<<"\t"<<*p;	
		}
		if(_cover[it->first].multimode>0)cout<<"\tMultimode="<<_cover[it->first].multimode;	
		if(PRINTCONSENSUP){
			ccl[ss.str()]=1;   ///ccl is a simple map with a string as key. Thus each neighbors set appears only once.
			if(ccl.size()>cclso){	///If the map increased in size it means that the last set of neighbors is a new one
				cclso=ccl.size(); ///update then the counter and add this set to the consensus partition.
				_consPart.insert(pair<int, string> ((it->second).neighbors.size(),ss.str()) );
			}
		}
		if(VERBOSE)if((it->second).neighbors.size()==1)cout<<"\tSingleton";
		cout<<endl;
	}
	cout<<"#EndCover"<<endl;
}

//void PartitionStats::getConsensusPartition(){
Partition PartitionStats::getConsensusPartition(){
	Partition _consensusPartition;
	if(!hasaCover()){
		cout<<"#W: Cover not yet defined. Setting up cover."<<endl;
		getCover();
	}
	if(!_consPart.size()>0){
		map<string,int> ccl;
		int cclso=0;
		for(map<string,sNeighborhood>::iterator it=_cover.begin();it!=_cover.end();it++){
			stringstream ss;
			for(sset::iterator p=(it->second).neighbors.begin();p!=(it->second).neighbors.end();p++){
				ss<<"\t"<<*p;	
			}
			ccl[ss.str()]=1;   ///ccl is a simple map with a string as key. Thus each neighbors set appears only once.
			if(ccl.size()>cclso){	///If the map increased in size it means that the last set of neighbors is a new one
				cclso=ccl.size(); ///update then the counter and add this set to the consensus partition.
				_consPart.insert(pair<int, string> ((it->second).neighbors.size(),ss.str()) );
			}
		}
	}
	_consensusPartition.xtrConsPart(_consPart);
	return _consensusPartition;
}

void PartitionStats::printConsensusPart(){
	int n=0;
	if(!QUIET)cout<<"#calculating Bell Number..."<<endl;
	//systemDate();
	//long int bn=4;
	//cout<<"#Bell_number ~ "<<BellN[_cover.size()]<<endl;
	//cout<<"#Bell_number ~ "<<BellN[bn]<<endl;
	//systemDate();
	if(!QUIET)cout<<"#BeginConsensusPartition"<<endl;
	for(multimap<int,string,greaterThan>::iterator cl=_consPart.begin();cl!=_consPart.end();cl++)
		cout<<cl->first<<"\tN"<<++n<<"\t"<<cl->second<<endl;
	if(!QUIET)cout<<"#EndConsensusPartition"<<endl;
}

void PartitionStats::setConsensus_Ad(){
	_Ad.clear();
	if(VERBOSE)cout<<"#PartitionStats::setConsensus_Ad() : _partitionl.size()="<<_partitionl.size()<<endl;
	ppvect::iterator p=_partitionl.begin();
	/* Ad_os() aims at displaying the adjacency matrix ordered in such a way that clusters locate in boxes around the diagonal.
	if(_partitionl.size()==1)
		_Ad=p->Ad_os();
	else
		_Ad=p->Ad();
	*/
	_Ad=p->Ad_os();
	for(p++ ;p!=_partitionl.end();p++){
		if(VERBOSE)cout<<"#PartitionStats::setConsensus_Ad() : "<<p->FileName();
		//graph g=p->Ad();
		graph g=p->Ad_os();
		_Ad+=g;
		if(VERBOSE)cout<<" ... OK"<<endl;
	}
	double z=1.0/_partitionl.size();
	_Ad=z*_Ad;
	if(_Ad.empty()){
		cout<<"ERROR: PartitionStats::setConsensus_Ad() : Failed to build Adjacency Matrix"<<endl;
		exit(1);
	}
}

void PartitionStats::get_Ad(){
	if(_Ad.empty()) setConsensus_Ad();
	cout<<"#BeginConsensusAdjacencyMatrix"<<endl;
	for(graph::iterator i=_Ad.begin();i!=_Ad.end();i++){
		string a = (i->first).first ;
		string b = (i->first).second ;
		strpair sp (a,b);
		cout<<a<<"\t"<<b<<"\t"<<_Ad[sp]<<endl;
	}
	cout<<"#EndConsensusAdjacencyMatrix"<<endl;
}

void  PartitionStats::get_FuzzyConsensusPartition(){
	if(VERBOSE) cout<<"#Setting up FuzzyConsensusPartition"<<endl;
	if(_Ad.empty()) setConsensus_Ad();
	vector< svect > fuzzyclusters;
	for(graph::iterator i=_Ad.begin();i!=_Ad.end();i++){
		string a = (i->first).first ;
		string b = (i->first).second ;
		strpair sp (a,b);
		if(_Ad[sp]>0 && a.compare(b)!=0){
			if(VERBOSE) cout<<"#Element "<<a<<" is in relation with "<<b<<endl;
			if(fuzzyclusters.size()==0){
				svect tes;
				tes.push_back(a);
				tes.push_back(b);
				fuzzyclusters.push_back(tes);
				if(VERBOSE) cout<<"#\tCreated initial fuzzy cluster"<<endl;
			}
			else {
				bool found=false;
				for(vector<svect >::iterator cl=fuzzyclusters.begin();cl!=fuzzyclusters.end();cl++){
					if( (find(cl->begin(),cl->end(),a)!=cl->end()) && (find(cl->begin(),cl->end(),b)==cl->end()) ){
						cl->push_back(b);
						found=true;
						if(VERBOSE) cout<<"#\tAdded to previous cluster of "<<a<<endl;
					}
					else if( (find(cl->begin(),cl->end(),b)!=cl->end()) && (find(cl->begin(),cl->end(),a)==cl->end()) ){
						cl->push_back(a);
						found=true;
						if(VERBOSE) cout<<"#\tAdded to previous cluster of "<<b<<endl;
					}
					else if( (find(cl->begin(),cl->end(),a)==cl->end()) && (find(cl->begin(),cl->end(),b)==cl->end()) ){
						found=true;
						if(VERBOSE) cout<<"#\tAlready present "<<b<<endl;
					}
					if(found)break;
				}
				if(!found){
					svect tes;
					tes.push_back(a);
					tes.push_back(b);
					fuzzyclusters.push_back(tes);
					if(VERBOSE) cout<<"#\tCreated a new fuzzy cluster"<<endl;
				}
			}
		}
	}
	containerLargerThan<svect> svectComparator;	
	sort(fuzzyclusters.begin(),fuzzyclusters.end(), svectComparator);
	for(vector<svect >::iterator cl=fuzzyclusters.begin();cl!=fuzzyclusters.end();cl++){
		//sort(cl->begin(),cl->end(),greaterThan);
		sort(cl->begin(),cl->end());
	}
	unique(fuzzyclusters.begin(),fuzzyclusters.end());
	cout<<"#BeginFuzzyConsensusPartition"<<endl;
	int n=0;
	for(int cl=0;cl<fuzzyclusters.size();cl++){
		cout<<fuzzyclusters[cl].size()<<"\tN"<<cl;
		for(svect::iterator sit=fuzzyclusters[cl].begin();sit!=fuzzyclusters[cl].end();sit++){
			cout<<"\t"<<*sit;
			n++;
		}
		cout<<endl;
	}
	cout<<"#EndFuzzyConsensusPartition"<<endl;
	//DETERMINE FUZZY CHARACTERISTIC FUNCTION
	vector< map<string, double > > characteristicFunction;
	characteristicFunction.resize(fuzzyclusters.size());
	//characteristicFunction.assign(fuzzyclusters.size(),sizeOf(map<string, double >));
	for(graph::iterator i=_Ad.begin();i!=_Ad.end();i++){
		string a = (i->first).first ;
		string b = (i->first).second ;
		svect tmpv;
		tmpv.assign(2,a);
		tmpv[1]=b;
		double mx=-1.0;
		for(svect::iterator ra=tmpv.begin();ra!=tmpv.end();ra++){
			for(int cl=0;cl<fuzzyclusters.size();cl++){
				for(svect::iterator sit=fuzzyclusters[cl].begin();sit!=fuzzyclusters[cl].end();sit++){
					if(ra->compare(*sit)==0)continue;
					strpair sp (*ra,*sit);
					if(_Ad[sp]>mx)mx=_Ad[sp];
					if(VERBOSE)cout<<"#\tmembership-a("<<*sit<<") of "<<*ra<<" into "<<cl<<" = "<<mx<<endl;
					sp=pair<string,string> (*sit,*ra);
					if(_Ad[sp]>mx)mx=_Ad[sp];
					if(VERBOSE)cout<<"#\tmembership-b("<<*sit<<") of "<<*ra<<" into "<<cl<<" = "<<mx<<endl;
				}
				if(VERBOSE)cout<<"#MEMBERSHIP of "<<*ra<<" into "<<cl<<" = "<<mx<<endl;
				characteristicFunction[cl].insert(pair<string,double> (*ra,mx) );
			}
		}
	}
	cout<<"#BeginFuzzyClusterMembership"<<endl;
	cout<<"#Clusters";
	for(int cl=0;cl<fuzzyclusters.size();cl++){
		cout<<"\tN"<<cl;
	}
	cout<<endl;
	//characteristicFunction vector<map<string::double>>
	for(map<string, double>::iterator el=characteristicFunction[0].begin();el!=characteristicFunction[0].end();el++){
		cout<<el->first; //first is the key.
		for(int cl=0;cl<fuzzyclusters.size();cl++){ //foreach cluster.
			cout<<"\t"<<characteristicFunction[cl][el->first]; //why not el->second?
		}
		cout<<endl;
	}	
	
	cout<<"#EndFuzzyClusterMembership"<<endl;

}

void PartitionStats::distancesPrintHeadComment(pmeasure measure, pmetricv metric, flagheader hd, bool usingREF){
	if(QUIET)return;
	string prefix="";
	switch(hd){
		case BEGIN:
			prefix="Begin";
			break;
		case END:
			prefix="End";
			break;
	}
	string entropy="";
	stringstream ss;
	switch(metric){
		//case entropy:
		case shannon:
			ss<<"Shannon";
			break;
		case cardinality:
			ss<<"Cardinality";
			break;
		case boltzmann:
			ss<<"Boltzman";
			break;
		case tsallis:
			ss<<"Tsallis extensivity= "<<_entensivity_degree;
			break;
		case renyi:
			ss<<"RenyiDistances extensivity= "<<_entensivity_degree;
			break;
		case jeffreyQnorm:
			ss<<"TarantolaJQNDistances extensivity= "<<_entensivity_degree;
			break;
	}
	entropy=ss.str();
	switch(measure){
		case conditionalEntropy:
			cout<<"#"<<prefix<<"ConditionalEntropy"<<": "<<entropy<<endl;
			if(!QUIET&&hd==BEGIN){
				if(!usingREF)cout<<"#PartitionA\tPartitionB\tH(A|B)\tH(B|A)\tVI-distance"<<endl;
				else {
					cout<<"#Partition(reference):"<<_partitionl[0].FileName()<<endl;
					cout<<"#PartitionB\tH(B|Ref)\tVI-distance"<<endl;
				}
			}
			break;
		case jointEntropy:
			cout<<"#"<<prefix<<"JointEntropy"<<": "<<entropy<<endl;
			if(!QUIET&&hd==BEGIN){
				if(!usingREF)
					cout<<"#PartitionA\tPartitionB\tH(A,B)"<<endl;
				else {
					cout<<"#Partition(reference):"<<_partitionl[0].FileName()<<endl;
					cout<<"#PartitionB\tH(B,Ref)"<<endl;
				}
			}
			break;
		case vmeasureArithmetic:
			cout<<"#"<<prefix<<"VmeasureArithmetic"<<": "<<entropy<<endl;
			if(!QUIET&&hd==BEGIN){
				if(!usingREF)
					cout<<"#PartitionA\tPartitionB\tVmeasure-Arithmetic"<<endl;
				else {
					cout<<"#Partition(reference):"<<_partitionl[0].FileName()<<endl;
					cout<<"#PartitionB\tVmeasure-Arithmetic"<<endl;
				}
			}
			break;
		case vmeasureGeometric: 
			cout<<"#"<<prefix<<"VmeasureGeometric"<<": "<<entropy<<endl;
			if(!QUIET&&hd==BEGIN){
				if(!usingREF)
					cout<<"#PartitionA\tPartitionB\tVmeasure-Geometric"<<endl;
				else {
					cout<<"#Partition(reference):"<<_partitionl[0].FileName()<<endl;
					cout<<"#PartitionB\tVmeasure-Geometric"<<endl;
				}
			}
			break;
		case vmeasureHarmonic:
			cout<<"#"<<prefix<<"VmeasureHarmonic"<<": "<<entropy<<endl;
			if(!QUIET&&hd==BEGIN){
				if(!usingREF)
					cout<<"#PartitionA\tPartitionB\tVmeasure-Harmonic"<<endl;
				else {
					cout<<"#Partition(reference):"<<_partitionl[0].FileName()<<endl;
					cout<<"#PartitionB\tVmeasure-Harmonic"<<endl;
				}
			}
			break;
		case symmetricPurity:
			cout<<"#"<<prefix<<"SymmetricPurityScores"<<endl;
			if(!QUIET&&hd==BEGIN){
				if(!usingREF)
					cout<<"#PartitionA\tPartitionB\tSymmetricPurityStrict\tSymmetricPurityLax"<<endl;
				else {
					cout<<"#Partition(reference):"<<_partitionl[0].FileName()<<endl;
					cout<<"#PartitionB\tSymmetricPurityStrict\tSymmetricPurityLax"<<endl;
				}
			}
			break;
	}
}

void PartitionStats::distancesPrintHeadComment(pmetricv metric, flagheader hd, bool usingREF){
	if(QUIET)return;
	string prefix="";
	switch(hd){
		case BEGIN:
			prefix="Begin";
			break;
		case END:
			prefix="End";
			break;
	}
	switch(metric){
		//case entropy:
		case shannon:
			cout<<"#"<<prefix<<"VIScores"<<endl;
			break;
		case cardinality:
			cout<<"#"<<prefix<<"EditScores"<<endl;
			break;
		case boltzmann:
			cout<<"#"<<prefix<<"BoltzmannDistances"<<endl;
			break;
		case tsallis:
			cout<<"#"<<prefix<<"TsallisDistances q= "<<_entensivity_degree<<endl;
			break;
		case renyi:
			cout<<"#"<<prefix<<"RenyiDistances q= "<<_entensivity_degree<<endl;
			break;
		case jeffreyQnorm:
			cout<<"#"<<prefix<<"TarantolaJQNDistances q= "<<_entensivity_degree<<endl;
			break;
	}
	if(hd==BEGIN&&!usingREF)cout<<"#Partition(target)\tPartition(ref)\tDist. target-intersection\tDist. ref-intersection\tDist. target-ref"<<endl;
}

void PartitionStats::distancesRef_Subsprojection(pmetricv metric){
	_DIST_SUBSPROJECT=true;
	if(!QUIET)cout<<"#Using DIST_SUBSPROJECT"<<endl;
	distancesRef(metric);
	_DIST_SUBSPROJECT=false;
}
void PartitionStats::distancesRef(pmetricv pm){
	Partition p=_partitionl[0];
	bool usingREF=true;
	distancesPrintHeadComment(pm,(flagheader) BEGIN, usingREF);
	if(!QUIET){
		cout<<"#Partition(reference):"<<p.FileName()<<endl;
		cout<<"#Partition(target)\ttarget-intersection\tref-intersection\ttarget-ref"<<endl;
	}
	for(int j=1;j<_partitionl.size();j++){
		if(_DIST_SUBSPROJECT) p=_partitionl[0];
		Partition pb=_partitionl[j];
		//if(p==pb)continue;
		Partition pi=p*pb;
		if(_DIST_SUBSPROJECT){
			if(VERBOSE)cout<<"#Projecting onto common subspace ("<<pi.sitems.size()<<") "<<endl;
			p.SubsProject(pi.sitems);
			pb.SubsProject(pi.sitems);
		}
#ifdef DEBUG
		cout<<"#_DIST_SUBSPROJECT="<<_DIST_SUBSPROJECT<<endl;
		cout<<"#"<<p.FileName()<<endl;
		p.printPartition();
		cout<<"#"<<pb.FileName()<<endl;
		pb.printPartition();
#endif
		switch(pm){
			//case entropy:
			case shannon:
				cout<<pb.FileName()<<"\t"<<VI(pb,pi)<<"\t"<<VI(pi,p)<<"\t"<<VI(pb,p)<<endl;
				break;
			case cardinality:
				cout<<pb.FileName()<<"\t"<<ES(pb,pi)<<"\t"<<ES(pi,p)<<"\t"<<ES(pb,p)<<endl;
				break;
			case boltzmann:
				cout<<pb.FileName()<<"\t"<<BK(pb,pi)<<"\t"<<BK(pi,p)<<"\t"<<BK(pb,p)<<endl;
				break;
			case tsallis:
				cout<<pb.FileName()<<"\t"<<TS(pb,pi)<<"\t"<<TS(pi,p)<<"\t"<<TS(pb,p)<<endl;
				break;
			case renyi:
				cout<<pb.FileName()<<"\t"<<RS(pb,pi)<<"\t"<<RS(pi,p)<<"\t"<<RS(pb,p)<<endl;
				break;
			case jeffreyQnorm:
				cout<<pb.FileName()<<"\t"<<TD(pb,pi)<<"\t"<<TD(pi,p)<<"\t"<<TD(pb,p)<<endl;
				break;
			default:
				cout<<"ERROR: PartitionStat::distancesRef() : unknown metric case"<<endl;
				exit(1);
		}
	}
	distancesPrintHeadComment(pm,(flagheader) END, usingREF);
}

void PartitionStats::distances_Subsprojection(pmetricv metric){
	_DIST_SUBSPROJECT=true;
	if(!QUIET)cout<<"#Using DIST_SUBSPROJECT"<<endl;
	distances(metric);
	_DIST_SUBSPROJECT=false;
}

void PartitionStats::distances(pmetricv pm){
	distancesPrintHeadComment(pm,(flagheader) BEGIN);
	for(int i=0;i<_partitionl.size()-1;i++){
		Partition pa=_partitionl[i];
		for(int j=i+1;j<_partitionl.size();j++){
			if(_DIST_SUBSPROJECT) pa=_partitionl[i];
			Partition pb=_partitionl[j];
			Partition pi=pa*pb;
			if(_DIST_SUBSPROJECT){
				if(VERBOSE)cout<<"#Projecting onto common subspace ("<<pi.sitems.size()<<") "<<endl;
				pa.SubsProject(pi.sitems);
				pb.SubsProject(pi.sitems);
			}
			switch(pm){
				//case entropy:
				case shannon:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<VI(pa,pi)<<"\t"<<VI(pi,pb)<<"\t"<<VI(pa,pb)<<endl;
					break;
				case cardinality:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<ES(pa,pi)<<"\t"<<ES(pi,pb)<<"\t"<<ES(pa,pb)<<endl;
					break;
				case boltzmann:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<BK(pa,pi)<<"\t"<<BK(pi,pb)<<"\t"<<BK(pa,pb)<<endl;
					break;
				case tsallis:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<TS(pa,pi)<<"\t"<<TS(pi,pb)<<"\t"<<TS(pa,pb)<<endl;
					break;
				case renyi:
					//cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<RS(pa,pi)<<"\t"<<RS(pi,pb)<<"\t"<<RS(pa,pb)<<endl;
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<RS(pa,pi)<<"\t"<<RS(pi,pb)<<"\t"<<RS(pa,pb)<<endl;
					break;
				case jeffreyQnorm:
					cout<<pa.FileName()<<"\t"<<pb.FileName()<<"\t"<<TD(pa,pi)<<"\t"<<TD(pi,pb)<<"\t"<<TD(pa,pb)<<endl;
					break;
				default:
					cout<<"ERROR: PartitionStat::distances() : unknown metric case"<<endl;
					exit(1);
			}
		}
	}	
	distancesPrintHeadComment(pm,(flagheader) END);
	_DIST_SUBSPROJECT=false;
}

void PartitionStats::getSplitsRef(splitmethod similarity=overlap){
	//cout.setf(ios::fixed,ios::floatfield);
	Partition* pa=&_partitionl[0];
	stringstream ss;
	stringstream sss;
	int grey;
	double value;
	cout<<"#BeginSplitAnalysis: Method=\t";
	switch(similarity){
		case cosine:	
			cout<<"Cosine Nomarlization Similarity [Leicht, E.A. et al., PRE, 73, 26120 (2006)]"<<endl;
			break;
		case split:
			cout<<"Split [Overlap fraction iif reference cluster is subset of target cluster; 0 otherwise]"<<endl;
			break;
		case overlap:
			cout<<"Overlap [Overlap fraction]"<<endl;
			break;
		default:
			cout<<"ERROR:  PartitionStats::getSplitRef : Unknown method."<<endl;
			exit(1);
			break;
	}
	cout<<"#Reference(columns): "<<pa->FileName()<<endl;
	for(int j=1;j<_partitionl.size();j++){
		cout<<"#Clusters\t";
		Partition* pb=&_partitionl[j];
		int width,height;
		width=pb->n_clusters(); 
		height=pa->n_clusters();
		if(similarity!=cosine){
			width++; //last column will contain the total for each cluster of pb
			height++; //last row will contain the total for each cluster of pb
		}
		ss<<"P2"<<endl;
		ss<<width<<" "<<height<<endl;
		ss<<_PGM_P2_GRAYSCALE_<<endl;
		for(smat::iterator cl=pb->clusters.begin();cl!=pb->clusters.end();cl++)
			cout<<*(cl->begin()+1)<<"\t";
		cout<<"Totals"<<endl;
		double atotal;
		int bidx;
		vector< double > btotal (pb->clusters.size(),0.0);
		vector< int > bsplit (pb->clusters.size(),0);
		for(smat::iterator cla=pa->clusters.begin();cla!=pa->clusters.end();cla++){
			cout<<*(cla->begin()+1)<<"\t";
			svect acl(cla->begin()+pa->cluster_offset(),cla->end());
			grey=0;
			value=atotal=0.0;
			bidx=-1;
			for(smat::iterator cl=pb->clusters.begin();cl!=pb->clusters.end();cl++){
				svect bcl(cl->begin()+pb->cluster_offset(),cl->end());
				svect incl=acl*bcl;
				value=1.0*incl.size();
				bidx++;
				switch(similarity){
					case cosine:
						value/=sqrt(1.0*acl.size()*bcl.size());
						break;
					case split:
						if(acl<=bcl) {
							value/=bcl.size() ;
						}
						else value=0.0;
						break;
					case overlap:
						//svect ucl=acl+bcl;
						value/=bcl.size();
						break;
				}
				grey=int( value*(_PGM_P2_GRAYSCALE_-1)+0.0 );
				printf("%4.2f\t",value);
				ss<<grey<<" ";
				//atotal+=value;//ROW TOTALS
				if(acl<=bcl){
					btotal[bidx]+=value; //COLUMN TOTALS
				}
				if(acl==bcl){
					atotal=1.0;//ROW TOTALS
				}
				if(incl.size()>0){
					bsplit[bidx]++;//NUMBER OF SPLITS FOR EACH COLUMN
				}
			}
			printf("%4.2f\n",atotal);//PRINT ROW TOTALS
			grey=int( atotal*(_PGM_P2_GRAYSCALE_-1)+0.0 );//PLOT ROW TOTALS
			ss<<grey;
			ss<<endl;
		}
		cout<<"Totals\t";
		for(int i=0;i<btotal.size();i++){
			printf("%4.2f\t",btotal[i]);//PRINT COLUMN TOTALS
			if(similarity!=cosine){
				grey=int( btotal[i]*(_PGM_P2_GRAYSCALE_-1)+0.0 );
				ss<<grey<<" ";//PLOT COLUMN TOTALS
			}
		}
		cout<<endl;
		if(similarity!=cosine){
			ss<<0<<endl;
		}
		cout<<"#NOverlaps\n#";
		for(smat::iterator cl=pb->clusters.begin();cl!=pb->clusters.end();cl++)
			cout<<*(cl->begin()+1)<<"\t";
		cout<<endl<<"#";
		for(int i=0;i<btotal.size();i++){
			printf("%5d\t",bsplit[i]);
		}
		cout<<endl;
	}
	cout<<"#EndSplitAnalysis"<<endl;
	cout<<"#BeginSplitImage"<<endl;
	cout<<ss.str();
	cout<<"#EndSplitImage"<<endl;
}

void PartitionStats::getPurityRef(){
	Partition* pa=&_partitionl[0];
	if(!QUIET)cout<<"#BeginPurityScores\n#Reference: "<<pa->FileName()<<endl;
	if(!QUIET)cout<<"#Partition(target)\tPurity Strict(raw)\tPurity Lax(raw)\tPurity Strict\tPurity Lax"<<endl;
	int ns,nl;
	double ps,pl;
	vector<double> pscores;
	for(int j=1;j<_partitionl.size();j++){
		Partition* pb=&_partitionl[j];
		pscores=pb->purityScore(pa);
		ps=pscores[0];
		pl=pscores[1];
		ns=(int)(ps*(double)pa->n_nonSingClusters());
		nl=(int)(pl*(double)pb->n_nonSingClusters());
		cout<<pb->FileName()<<"\t"<<ns<<"\t"<<nl<<"\t"<<ps<<"\t"<<pl<<endl;
	}
	if(!QUIET)cout<<"#EndPurityScores"<<endl;
}

void PartitionStats::getPurity(){
	if(!QUIET)cout<<"#BeginPurityScores"<<endl;
	if(!QUIET)cout<<"#Partition(target)\tPartition(ref)\tPurity Strict(raw)\tPurity Lax(raw)\tPurity Strict\tPurity Lax"<<endl;
	int ns,nl;
	double ps,pl;
	vector<double> pscores;
	for(int i=0;i<_partitionl.size()-1;i++){
		Partition* pa=&_partitionl[i];
		for(int j=i+1;j<_partitionl.size();j++){
			Partition* pb=&_partitionl[j];
			pscores=pa->purityScore(pb);
			ps=pscores[0];
			pl=pscores[1];
			ns=(int)(ps*(double)pb->n_nonSingClusters());
			nl=(int)(pl*(double)pa->n_nonSingClusters());
			cout<<pa->FileName()<<"\t"<<pb->FileName()<<"\t"<<ns<<"\t"<<nl<<"\t"<<ps<<"\t"<<pl<<endl;
			pscores=pb->purityScore(pa);
			ps=pscores[0];
			pl=pscores[1];
			ns=(int)(ps*(double)pa->n_nonSingClusters());
			nl=(int)(pl*(double)pb->n_nonSingClusters());
			cout<<pb->FileName()<<"\t"<<pa->FileName()<<"\t"<<ns<<"\t"<<nl<<"\t"<<ps<<"\t"<<pl<<endl;
		}
	}
	if(!QUIET)cout<<"#EndPurityScores"<<endl;
}

void PartitionStats::printHasseDiagram(){
	long int level,nextlevel;
	map<long int, ppvect >::reverse_iterator rit=_hasseNodes.rbegin();
	long int maxlevel=rit->first;
	cout<<"#BeginHasseDiagram"<<endl;
	cout<<"#Partition1  >  Partition2\tedge\tPart1-#clusters\tPart2-#nclusters"<<endl;
	map<long int, ppvect >::iterator h1=_hasseNodes.begin();
	level=h1->first;
	for( ;level!=maxlevel; h1++, level=h1->first){
		map<long int, ppvect >::iterator h2=h1;
		h2++;
		nextlevel=h2->first;
		if(VERBOSE)cout<<"#level="<<level<<" nextlevel="<<nextlevel<<" ("<<maxlevel<<") Checking intersection:"<<endl;
		for(ppvect::iterator pa=_hasseNodes[h1->first].begin();pa!=_hasseNodes[h1->first].end();pa++)
			for(ppvect::iterator pb=_hasseNodes[h2->first].begin();pb!=_hasseNodes[h2->first].end();pb++){
				Partition q=(*pa)*(*pb);
				if(VERBOSE)cout<<"#pa="<<(*pa).FileName()<<" /\\ "<<"pb="<<(*pb).FileName()<<endl;
					if(q==*pb){
						cout<<(*pa).FileName()<<"\t"<<(*pb).FileName()<<"\t"<<"1"<<"\t"<<(*pa).n_clusters()<<"\t"<<(*pb).n_clusters()<<endl;
					}
			}
	}
	cout<<"#EndHasseDiagram"<<endl;
}

void PartitionStats::printHasseNodes(){
	cout<<"#BeginHasseNodes (hierarchy)"<<endl;
	for(map<long int, ppvect >::iterator h=_hasseNodes.begin();h!=_hasseNodes.end();h++){
		cout<<h->first<<"\t"<<(h->second).size();
		for(ppvect::iterator pa=_hasseNodes[h->first].begin();pa!=_hasseNodes[h->first].end();pa++)
			cout<<"\t"<<pa->FileName();
		cout<<endl;
	}
	cout<<"#EndHasseNodes"<<endl;
}

#endif //END  _CLASS_PARTITIONSTATS
