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


#ifndef _CLASS_MULTIPLESEQALIGN
#define _CLASS_MULTIPLESEQALIGN 1

#include "MultipleSeqAlign.h"

MultipleSeqAlign::MultipleSeqAlign(){
        _msaf=NULL;
        _len=_nseq=0;
}

MultipleSeqAlign::MultipleSeqAlign(char* msaf)
{
        _msaf=msaf;
        _len=_nseq=0;
        ifstream _is(_msaf);
        if(!_is){
                cout<<"ERROR: MultipleSeqAlign : File not found: "<<_msaf<<endl;
                exit(1);
        }
        if(!QUIET) cout<<"#Reading MSA "<<msaf<<endl;
        string str; //String read from file
        string seq=""; //String containing sequence read from file
        Sequence Seq; //Actual sequence structure
        while(_is>>str){
                string fchr=str.substr(0,1);
                if(fchr.compare("#")==0)continue;
                if(VERBOSE)cout<<"#read: "<<str<<" First character="<<fchr<<endl;
                if(strcmp(fchr.c_str(),">")==0){
                        if(VERBOSE)Seq.print();
                        if(VERBOSE)cout<<"#New string: "<<str<<endl;
                        if(!seq.compare("")==0){
                                Seq.setSeq(seq);
                                if(_len==0)_len=Seq.length();
                                if(!Seq.length()==_len){
                                        cout<<"ERROR: Sequence length doesn't fit MSA length ("<<_len<<") : "<<Seq.name()<<" ("<<Seq.length()<<")"<<endl;
                                        exit(1);
                                }
                                _Seqlist.push_back(Seq);
                        }
                        Seq.setName(str);
                        seq="";
                        continue;
                }
                seq+=str;
        }
        Seq.setSeq(seq);        ///The last sequence isn't followed by a line starting with ">"
        _Seqlist.push_back(Seq);
        _nseq=_Seqlist.size();
        if(!QUIET) cout<<"#Found "<<_nseq<<" sequences of length "<<_len<<". Last read: "<<(_Seqlist.end()-1)->name()<<endl;
}

MultipleSeqAlign::MultipleSeqAlign(MSA& msa){
	if(! msa.size()>0){
		cout<<"ERROR: Cannot build MSA"<<endl;
		exit(1);
	}
	_msaf=NULL;
	_Seqlist=msa;
	_len=_Seqlist[0].length();
	_nseq=_Seqlist.size();
}

void MultipleSeqAlign::addSeq(Sequence Seq){
	addSeq(&Seq);
}

void MultipleSeqAlign::addSeq(Sequence* Seq){
        if(_len>0)
		if(_len!=Seq->length()){
                	cout<<"Sequence length "<<Seq->name()<<" ("<<Seq->length()<<") doesn't fit MSA's length ("<<_len<<")"<<endl;
                	exit(1);
        	}
        else {
                _len=Seq->length();
                _nseq=0;
        }
        _Seqlist.push_back(*Seq);
        _nseq++;
}

MultipleSeqAlign  MultipleSeqAlign::xtractPositions(vector<int>* positions){
	MultipleSeqAlign nmsa;
	nmsa.setName("NAA");
	if(positions->empty())return nmsa;
	string sq;
	string name;
	for(MSA::iterator st=_Seqlist.begin();st!=_Seqlist.end();st++){
		name=st->name();
		sq=st->xtractPositions(positions);
		Sequence seq(name,sq);
		nmsa.addSeq(&seq);
	}
	if(!nmsa.empty())nmsa.setName(getName());
	return nmsa;
}

MultipleSeqAlign  MultipleSeqAlign::xtractSequences(svect* seqnames, bool equal){
	MultipleSeqAlign nmsa;
	string name;
	nmsa.setName("Extracted sequences");
	bool found,cond;
	for(MSA::iterator st=_Seqlist.begin();st!=_Seqlist.end();st++){
		name=st->name();
		found=false;
		for(svect::iterator nit=seqnames->begin();nit!=seqnames->end();nit++){
			cond=(name.compare(*nit)==0||name.substr(1,name.length()-1).compare(*nit)==0); //Name may start by >
			if(cond){
				found=true;
				if(equal)nmsa.addSeq(*st); //If equal, add sequences found in list; else
			}
		}
		if(!(equal||found))nmsa.addSeq(*st); //... add sequences NOT in list
	}
	if(nmsa.empty()&&!QUIET)cout<<"WARNING : NO sequence extracted"<<endl;
	return nmsa;
}
///For each sequence of MSA-B, choose among top 10 highest ID hits against MSA-A, the top one that hasn't been chosen yet.
MultipleSeqAlign  MultipleSeqAlign::xtractSequencesHighestId(MultipleSeqAlign* msab, int cullsize, vector<int>* positions){
	MultipleSeqAlign nmsa;
	string seqstr;
	nmsa.setName("Extracted sequences");
	bool cond;
	double sid;
	double minId=0.0;
	double thrId=0.0;
	int oldsz=0;
	int seq_oldsz=0;
	int TOPLISTSIZE=10;
	map<string,Sequence> non_redundant_list; //Non-redundant list of top-id sequences found so far
	multimap<double,Sequence,greaterThan> topid; //Same as non_redundant_list but ordered by sequence ID starting from highest.
	for(MSA::iterator stb=msab->beginSeq();stb!=msab->endSeq();stb++){
		map<string,Sequence> seq_non_redundant_list;
		multimap<double,Sequence,greaterThan> st_topid;
		for(MSA::iterator st=_Seqlist.begin();st!=_Seqlist.end();st++){
			seqstr=st->sequence();
			sid=st->id(*stb,positions);
			cond=(sid>minId);
			if(cond||oldsz<TOPLISTSIZE){ //Cull top TOPLISTSIZE sequences from MSA-A
				if(seq_oldsz==TOPLISTSIZE){
					seq_non_redundant_list.erase(seq_non_redundant_list.rend()->first); //remove last item
				}
				seq_non_redundant_list[seqstr]=*st;
				if(seq_oldsz<seq_non_redundant_list.size()){ //Update lists if sequence hasn't been added before
					seq_oldsz=seq_non_redundant_list.size(); //update number of top-id sequences found so far. 
					st_topid.insert(pair<double,Sequence>(sid,*st)); //update top id list . We could use this map as well for updating oldsz
					minId=st_topid.rend()->first; //Update lowest top id value
				}
			}
		}
		for(multimap<double,Sequence,greaterThan>::iterator it=st_topid.begin();it!=st_topid.end();it++){
			non_redundant_list[(it->second).sequence()]=it->second; //Try to add sequence to non-redundant list
			if(oldsz<non_redundant_list.size()){ //... and if it hasn't been added before
				oldsz=non_redundant_list.size();
				topid.insert(pair<double,Sequence>(it->first,it->second));
				break; //One sequence is enough
			}
		}
	}
	int i=1;
	for(multimap<double,Sequence,greaterThan>::iterator it=topid.begin();it!=topid.end();it++,i++){
		nmsa.addSeq(it->second);  //... add sequence to new MSA.
		if(cullsize>0&&i==cullsize) break; //Stop searching for sequences if we found already cullsize
	}
	if(nmsa.empty()&&!QUIET)cout<<"#WARNING : NO sequence extracted"<<endl;
	return nmsa;
}

MultipleSeqAlign  MultipleSeqAlign::xtractSequencesById(MultipleSeqAlign* msab, double minId, double maxId, bool equal,vector<int>* positions){
	MultipleSeqAlign nmsa;
	string seqstr;
	nmsa.setName("Extracted sequences");
	bool found,cond;
	double sid;
	int oldsz;
	map<string,Sequence> non_redundant_list;
	for(MSA::iterator st=_Seqlist.begin();st!=_Seqlist.end();st++){
		oldsz=non_redundant_list.size();
		seqstr=st->sequence();
		found=false;
		for(MSA::iterator stb=msab->beginSeq();stb!=msab->endSeq();stb++){
			sid=st->id(*stb,positions);
			cond=(sid>minId&&sid<=maxId);
			if(cond){
				found=true;
				if(equal){
					non_redundant_list[seqstr]=*st; //If equal, add sequences found in list; else
					break; //One hit is enough
				}
			}
		}
		if(!(equal||found)){
			non_redundant_list[seqstr]=*st; //... add sequences NOT in list
				//One hit is enough NOT ENOUGH to be chosen
		}
		if(oldsz<non_redundant_list.size()){ //Add sequence if it hasn't been added before
			nmsa.addSeq(*st);
			oldsz=non_redundant_list.size();
		}
	}
	if(nmsa.empty()&&!QUIET)cout<<"#WARNING : NO sequence extracted"<<endl;
	return nmsa;
}

double MultipleSeqAlign::SeqId(int Seqn, int Seqm, vector<int>* positions)
{
        if(Seqn<0||Seqm<0||Seqn>_nseq||Seqm>_nseq){
                cout<<"Sequence indexes out of bound : "<<Seqn<<" , "<<Seqm<<endl;
                exit(1);
        }
        return _Seqlist[Seqn].id( _Seqlist[Seqm] , positions);
}

double MultipleSeqAlign::averageId(vector<int>* positions)
{
        if(!QUIET)cout<<"#Starting averageId calculation"<<endl;
        double avid=0.0;
        if(VERBOSE){
		cout<<"#_Seqlist.size()="<<_Seqlist.size();
		if(positions!=NULL)cout<<" positions->empty()"<<positions->empty()<<endl;
		else cout<<" positions is NULL"<<endl;
	}
        for(MSA::iterator ita=_Seqlist.begin();ita!=_Seqlist.end()-1;ita++){
                for(MSA::iterator itb=ita+1;itb!=_Seqlist.end();itb++){
                        if(VERBOSE)cout<<"#calculating pair-wise id..."<<endl;
                        avid+=(*ita).id(*itb,positions);
                        if(VERBOSE){
                                cout<<"# "<<(*ita).name()<<" - "<<(*itb).name()<<" id="<<(*ita).id(*itb,positions)<<endl;
                                (*ita).print();
                                (*itb).print();
                                for(int i=0;i<_len;i++)
                                        if((i+1)%10>0)
                                                cout<<(i+1)%10;
                                        else
                                                cout<<" ";
                                cout<<endl;
                        }
                }
	}
        return avid/(double)(_nseq*(_nseq-1)*0.5);
}

void MultipleSeqAlign::printPairwiseIds(MultipleSeqAlign& msa, vector<int>* positions)
{
        for(MSA::iterator ita=_Seqlist.begin();ita!=_Seqlist.end();ita++)
                for(MSA::iterator itb=msa.beginSeq();itb!=msa.endSeq();itb++){
                        string na=(*ita).name();
                        if(na.substr(0,1).compare(">")==0)na=na.substr(1,na.length()-1);
                        string nb=(*itb).name();
                        if(nb.substr(0,1).compare(">")==0)nb=nb.substr(1,nb.length()-1);
                        cout<<na<<"\t"<<nb<<"\t"<<SeqId(*ita,*itb,positions)<<endl;
                }
}

void MultipleSeqAlign::printPairwiseIds(vector<int>* positions)
{
        for(MSA::iterator ita=_Seqlist.begin();ita!=_Seqlist.end()-1;ita++)
                for(MSA::iterator itb=ita+1;itb!=_Seqlist.end();itb++){
                        string na=(*ita).name();
                        if(na.substr(0,1).compare(">")==0)na=na.substr(1,na.length()-1);
                        string nb=(*itb).name();
                        if(nb.substr(0,1).compare(">")==0)nb=nb.substr(1,nb.length()-1);
                        cout<<na<<"\t"<<nb<<"\t"<<SeqId(*ita,*itb,positions)<<endl;
                        //cout<<(*ita).name()<<"\t"<<(*itb).name()<<"\t"<<SeqId(*ita,*itb)<<endl;
                }
}

void MultipleSeqAlign::printAveragePairwiseIds(MultipleSeqAlign& msa, double thr, vector<int>* positions)
{
	cout<<"#Id_threshold= "<<thr<<endl;
        int onn=0;
        //double thr=50.0;
        for(MSA::iterator ita=_Seqlist.begin();ita!=_Seqlist.end();ita++){
                string na=(*ita).name();
                if(na.substr(0,1).compare(">")==0)na=na.substr(1,na.length()-1);
                double id,id2,v,mdid;
                id=id2=v=mdid=0.0;
                double idmin=101;
                double idmax=-1;
                int nn=0;
                int nthr=0;
                for(MSA::iterator itb=msa.beginSeq();itb!=msa.endSeq();itb++){
                        if(ita==itb)continue;
                        v=SeqId(*ita,*itb,positions);
                        v=v<=1?v*100:v;
                        if(v>thr)nthr++;
                        id+=v;
                        id2+=v*v;
                        idmin=v<idmin?v:idmin;
                        idmax=v>idmax?v:idmax;
                        nn++;
                }
                ///Print AverageID  StdDev. Variance Min Max #edges>50%Id %edges>50 total#edges
                cout<<na<<"\t"<<id/(double)nn<<"\t"<<sqrt(id2/nn-(id/nn)*(id/nn))<<"\t"<<(id2/nn-(id/nn)*(id/nn))<<"\t"<<idmin<<"\t"<<idmax<<"\t"<<nthr<<"\t"<<(double)nthr/nn*100.0<<"\t"<<nn<<endl;
                if(VERBOSE&&!onn==0&&!onn==nn)
                        cout<<"#WARNING: nodes with different number of edges"<<endl;
                onn=nn;
        }
}

void MultipleSeqAlign::printAveragePairwiseIds(double thr, vector<int>* positions)
{
	cout<<"#Id_threshold= "<<thr<<endl;
        int onn=0;
        //double thr=50.0;
        for(MSA::iterator ita=_Seqlist.begin();ita!=_Seqlist.end()-0;ita++){
                string na=(*ita).name();
                if(na.substr(0,1).compare(">")==0)na=na.substr(1,na.length()-1);
                double id,id2,v,mdid;
                id=id2=v=mdid=0.0;
                double idmin=101;
                double idmax=-1;
                int nn=0;
                int nthr=0;
                for(MSA::iterator itb=_Seqlist.begin();itb!=_Seqlist.end();itb++){
                        if(ita==itb)continue;
                        v=SeqId(*ita,*itb,positions);
                        v=v<=1?v*100:v;
                        if(v>thr)nthr++;
                        id+=v;
                        id2+=v*v;
                        idmin=v<idmin?v:idmin;
                        idmax=v>idmax?v:idmax;
                        nn++;
                }
                ///Print AverageID  StdDev. Variance Min Max #edges>50%Id %edges>50 total#edges
                cout<<na<<"\t"<<id/(double)nn<<"\t"<<sqrt(id2/nn-(id/nn)*(id/nn))<<"\t"<<(id2/nn-(id/nn)*(id/nn))<<"\t"<<idmin<<"\t"<<idmax<<"\t"<<nthr<<"\t"<<(double)nthr/nn*100.0<<"\t"<<nn<<endl;
                if(VERBOSE&&!onn==0&&!onn==nn)
                        cout<<"#WARNING: nodes with different number of edges"<<endl;
                onn=nn;
        }
}

void MultipleSeqAlign::print()
{
        for(MSA::iterator ita=_Seqlist.begin();ita!=_Seqlist.end()-0;ita++)
                (*ita).print();
}

MultipleSeqAlign MultipleSeqAlign::genRedundantMSA(int n){
/*
	unsigned int seed=time(NULL);
	genRedundantMSA(n,seed);
}
MultipleSeqAlign MultipleSeqAlign::genRedundantMSA(int n, unsigned int seed){
*/
	if(_Seqlist.size()==0){ 
		cout<<"ERROR: No MSA defined"<<endl;
		exit(1);
	}
	vector< int> randseq=ruleta.spinWheeli(n,_Seqlist.size());
	MSA redlist=_Seqlist;
	map<int, int> dg;
	stringstream redstr;
	//unsigned long int tr=0;
	long double tr=0.0;
#ifdef DEBUG
	cout<<"#DEBUG: picking up random seq..."<<endl;
#endif
	for(int i=0;i<randseq.size();i++){
#ifdef DEBUG
		cout<<"#DEBUG: Seq# "<<randseq[i]<<endl;
#endif
		Sequence rs=_Seqlist[randseq[i]];
		if(redstr.str().size()>0)redstr<<"\n"<<randseq[i];
		else redstr<<randseq[i];
		tr+=pow((long double)randseq[i],(long double)2.0);
		if(dg.find(randseq[i])==dg.end())dg.insert(pair<int,int>(randseq[i],0));
		stringstream ss;
		ss<<_Seqlist[randseq[i]].name()<<"-RED-"<<++dg[randseq[i]];
		rs.setName(ss.str());
		redlist.push_back(rs);
	}
#ifdef DEBUG
	cout<<"#DEBUG: done"<<endl;
#endif
	MultipleSeqAlign red(redlist);
	red.setName(redstr.str());
	tr=sqrt(tr);
	stringstream sstr;
	sstr<<"R"<<n<<"-"<<tr;
	red.setFileName((char*)sstr.str().c_str());
	return red;
}

void MultipleSeqAlign::NsamplesRedMSA(int nsamples=1, int n=1){
	unsigned int seed=time(NULL);
	NsamplesRedMSA(seed, nsamples, n);
}
void MultipleSeqAlign::NsamplesRedMSA(unsigned int seed, int nsamples=1, int n=1){
	ruleta.setSeed(seed);
	if(!QUIET){
		 cout<<"#Generating "<<nsamples<<" samples of MSAs with "<<n<<" sequences redundancy"<<endl;
		cout<<"#Random seed: "<<seed<<endl;
		cout<<"#Samples: ";
	}
	for(int i=0;i<nsamples;i++){
		if(!QUIET) cout<<i+1<<"-";
		MultipleSeqAlign redmsa=genRedundantMSA(n);
		stringstream redstr;
		char* fn=redmsa.getFileName();
		if(!QUIET) cout<<fn<<" ";
		redstr<<fn<<"/"<<fn<<".fasta";
#ifdef DEBUG
		cout<<"#DEBUG: Done first sample...writing..."<<endl;
		cout<<"#DEBUG: path : "<<redstr.str()<<endl;
#endif
		DIR* pdir=opendir(fn);
		if(!pdir){
			mode_t perm=0744;
			int ci=mkdir(fn,perm);
			if(ci!=0){
				cout<<"ERROR: could not create directory "<<fn<<endl;
				exit(1);
			}
		}
		else if(!QUIET) cout<<"\n#WARNING: directory "<<fn<<" already exists. \n#Samples: ";
			/*
		errno=0;
		struct dirent* pent;
		cout<<"#DEBUG: checking if exists file "<<redstr.str().c_str()<<endl;
		while((pent=readdir(pdir))){
			cout<<"#WARNING: directory "<<fn<<" already contains file "<<pent->d_name<<endl;
			//cout<<"#WARNING: directory "<<fn<<" already contains file "<<endl;
			if(memcmp(pent->d_name,redstr.str().c_str(),sizeof(pent->d_name))==0){
				cout<<"ERROR: File "<<redstr<<" already exists!"<<endl;
				exit(1);
			}
		}
		if(errno){
			cout<<"ERROR: failure while reading directory "<<fn<<endl;
			exit(1);
		}
			*/
#ifdef DEBUG
		cout<<"#DEBUG: writing .fasta file "<<redstr.str().c_str()<<endl;
#endif
		ofstream of(redstr.str().c_str());
		of<<redmsa;
		of.close();
		closedir(pdir);
#ifdef DEBUG
		getchar();
#endif
	}
	if(!QUIET) cout<<endl;
	ruleta.writeSeed();
}

ostream& operator<<(ostream& os, MultipleSeqAlign& msa){
	for(MSA::iterator seq=msa._Seqlist.begin();seq!=msa._Seqlist.end();seq++)
		os<<*seq<<endl;
	return os;
}

#endif //END _CLASS_MULTIPLESEQALIGN
