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
	if(!QUIET)printPIDNormalization();
	
}

MultipleSeqAlign::MultipleSeqAlign(char* msaf)
{
        _msaf=msaf;
        _len=_nseq=0;
	if(!QUIET)printPIDNormalization();
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
                if(fchr.compare("#")==0){
				getline(_is,str);
				continue;
		}
                //Found beginning of new sequence: push old one to the MSA and set name of the current new one
                if(strcmp(fchr.c_str(),">")==0){
                        if(!seq.compare("")==0){ //If we have the actual sequence...
                                Seq.setAlignedSeq(seq); //set it...
                                if(_len==0)_len=Seq.alignmentLength();///and check whether it may belong or not to the same MSA
                                if(!Seq.alignmentLength()==_len){
                                        cout<<"ERROR: Sequence length doesn't fit MSA length ("<<_len<<") : "<<Seq.name()<<" ("<<Seq.alignmentLength()<<")"<<endl;
                                        exit(1);
                                }
                                _Seqlist.push_back(Seq); ///If the length fits, add the previous sequence Seq to the MSA _Seqlist
                        	if(DEBUG)Seq.printAlignment();
                        }
                        if(DEBUG)cout<<"#New string: "<<str<<endl;
                        Seq.setName(str);//Set name of the current new sequence found
                        seq="";
                        continue;
                }
                if(DEBUG)cout<<"#read: "<<str<<" First character="<<fchr<<endl;
                seq+=str;
        }
        Seq.setAlignedSeq(seq);        ///The last sequence isn't followed by a line starting with ">"
	if(DEBUG)Seq.printAlignment();
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
	_len=_Seqlist[0].alignmentLength();
	_nseq=_Seqlist.size();
}

void MultipleSeqAlign::addSeq(Sequence Seq){
	addSeq(&Seq);
}

void MultipleSeqAlign::addSeq(Sequence* Seq){
        if(_len>0)
		if(_len!=Seq->alignmentLength()){
                	cout<<"Sequence length "<<Seq->name()<<" ("<<Seq->alignmentLength()<<") doesn't fit MSA's length ("<<_len<<")"<<endl;
                	exit(1);
        	}
        else {
                _len=Seq->alignmentLength();
                _nseq=0;
        }
        _Seqlist.push_back(*Seq);
        _nseq++;
}

//Extract sequences specified by name. IF bool equal=false, drop them instead.
MultipleSeqAlign  MultipleSeqAlign::xtractSequences(svect* seqnames, bool equal){
	MultipleSeqAlign nmsa;
	string name;
	nmsa.setName("Extracted sequences");
	bool found,cond;
	for(MSA::iterator st=_Seqlist.begin();st!=_Seqlist.end();st++){
		name=st->name();
		found=false;
		if(VERBOSE)cout<<name<<"? :";
		for(svect::iterator nit=seqnames->begin();nit!=seqnames->end();nit++){
			cond=(name.compare(*nit)==0||name.substr(1,name.length()-1).compare(*nit)==0); //Name may start by >
			if(cond){
				found=true;
				if(equal){
					nmsa.addSeq(*st); //If equal, add sequences found to list; else
					if(VERBOSE)cout<<" (Equal="<<equal<<"/found="<<found<<") Cull :"<<*nit;
				}else if(VERBOSE)cout<<" (Equal="<<equal<<"/found="<<found<<")  Drop :"<<*nit;
			}
		}
		if(!(equal||found)){
			nmsa.addSeq(*st); //... add sequences NOT yet in list
			if(VERBOSE)cout<<" Keeping ";
		}else if(VERBOSE)cout<<" Ignoring :";
		if(VERBOSE)cout<<endl;
	}
	if(nmsa.empty()&&!QUIET)cout<<"WARNING : NO sequence extracted"<<endl;
	if(!QUIET)cout<<"#Extracted "<<nmsa.getNumberOfSeq()<<" sequences"<<endl;
	return nmsa;
}

Partition MultipleSeqAlign::getPartition(char* msaf,MSAformat fmt){
	int offset=2;
	bool dosort=true;
	//char* partf;
	//sprintf(partf,"%s lst",msaf);
	smat pclustersl;
	svect cluster;
        ifstream _is(msaf);
        if(!_is){
                cout<<"ERROR: MultipleSeqAlign::getPartition : File not found: "<<msaf<<endl;
                exit(1);
        }
        if(!QUIET) cout<<"#Reading MSA "<<msaf<<endl;
        string str; //String read from file
	string cln="";
	long int nelements=0;
        while(_is>>str){
                string fchr=str.substr(0,1);
                if(fchr.compare("#")==0){
				getline(_is,str);
				continue;
		}
                if(fchr.compare("%")==0||\
		   fchr.compare("=")==0\
		  ){
			if(DEBUG)cout<<"#FOUND GROUP: "<<str<<" with name ";
			if(cluster.size()>0){
				stringstream ss;
				ss<<nelements;
				cluster[0]=ss.str();
				pclustersl.push_back(cluster);
				cluster.clear();
				nelements=0;
			}
			cln=str.substr(1,str.length()-1);
			while(cln.substr(0,1).compare("=")==0||cln.substr(0,1).compare(" ")==0){
				cln=cln.substr(1,str.length()-1);
			}
			cluster.push_back("0"); //Number of elements
			cluster.push_back(cln); //Name of cluster
			if(DEBUG){
				cout<<cln<<"\n#Current number of clusters= "<<pclustersl.size()<<endl;
			}
			getline(_is,str);
			continue;
		}
                //Found beginning of new sequence: push old one to the MSA and set name of the current new one
                if(strcmp(fchr.c_str(),">")==0){
                        if(DEBUG)cout<<"#New element: "<<str<<endl;
			nelements++;
                        cluster.push_back(str.substr(1,str.length()-1));
                        continue;
                }
        }
	stringstream ss;
	ss<<nelements;
	cluster[0]=ss.str();
	pclustersl.push_back(cluster);
	cluster.clear();
	//return Partition(&pclustersl,offset,dosort,partf);
	return Partition(&pclustersl,offset,dosort,msaf);
}

void MultipleSeqAlign::printWithClusterLabels(Partition* part, MSAformat fmt){
	string seqn_suffix="";
	int skipped=0;
	if(fmt==msaFmtNULL)fmt=GSIM;
	stringstream clsizes;
	if(!QUIET)cout<<"#BeginMSAwClusters"<<endl;
	for(smat::iterator cl=part->clusters.begin();cl!=part->clusters.end();cl++){
		switch(fmt){
			case FASTA3:
			case GDE3:
				if(cl->size()-part->cluster_offset()<3){
					//if(!QUIET)cout<<"#Skipping cluster with <3 elements: "<<part->getClusterName(*cl)<<endl;
					skipped++;
					continue;
				}
			case FASTA2:
			case GDE2:
				if(cl->size()-part->cluster_offset()<2){
					skipped++;
					continue;
				}
			case FASTA:
			case GDE:
				cout<<"==cluster_"<<part->getClusterName(*cl)<<" =="<<endl;
				break;
			//#####################
			case GSIM3:
				if(cl->size()-part->cluster_offset()<3){
					skipped++;
					continue;
				}
			case GSIM2:
				if(cl->size()-part->cluster_offset()<2){
					skipped++;
					continue;
				}
			case GSIM:
				seqn_suffix=part->getClusterName(*cl);
				break;
			//#####################
			case SPEER3:
				if(cl->size()-part->cluster_offset()<3){
					skipped++;
					continue;
				}
			case SPEER2:
				if(cl->size()-part->cluster_offset()<2){
					skipped++;
					continue;
				}
			case SPEER:
				clsizes<<cl->size()-part->cluster_offset()<<" ";
				break;
			//#####################
		}
		for(svect::iterator it=cl->begin()+part->cluster_offset();it!=cl->end();it++){
			for(MSA::iterator st=_Seqlist.begin();st!=_Seqlist.end();st++){
				string clit,seqn;
				clit=*it;
				seqn=st->name();
				if(clit.compare(seqn)==0||clit.compare(seqn.substr(1,seqn.length()-1))==0\
				   ||clit.substr(1,clit.length()-1).compare(seqn)==0\
				   ||clit.substr(1,clit.length()-1).compare(seqn.substr(1,seqn.length()-1))==0\
					)st->printAlignment(fmt,seqn_suffix);
			}	
		}
	}
	if(fmt==SPEER||fmt==SPEER3||fmt==SPEER2){
		if(!QUIET)cout<<"#ClusterSizes"<<endl;
		cout<<clsizes.str()<<endl;
	}
	if(!QUIET&&(fmt==FASTA2||fmt==FASTA3||fmt==GDE2||fmt==GDE3||fmt==SPEER2||fmt==SPEER3))cout<<"#Clusters shown: "<<part->clusters.size()-skipped<<endl;
	if(!QUIET)cout<<"#EndMSAwClusters"<<endl;
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

MultipleSeqAlign  MultipleSeqAlign::dropClones(MultipleSeqAlign* msab){
	MultipleSeqAlign nmsa;
	string msaname="Noclone"+msab->getName();
	nmsa.setName(msaname);
	map<string,Sequence> non_redundant_list; //Non-redundant list of top-id sequences found so far
	int oldsz=0;
	bool keep;
	for(MSA::iterator st=beginSeq();st!=endSeq();st++){
		keep=true;
		for(MSA::iterator stb=msab->beginSeq();stb!=msab->endSeq();stb++){
			if(st->alignedSequence().compare(stb->alignedSequence())!=0)continue;//If st is not a clone of stb, check against next stb
			keep=false; //If st is a clone of stb ...
			break;
		}
		if(keep){ //... do not keep st and check the next st
			non_redundant_list[st->alignedSequence()]=*st;
			if(oldsz<non_redundant_list.size()){ //If we didn't cull st before, ...
				nmsa.addSeq(*st); //... add this sequence.
				oldsz=non_redundant_list.size();
			}
		}
	}
	if(nmsa.empty()&&!QUIET)cout<<"#WARNING : NO sequence extracted"<<endl;
	return nmsa;
}

MultipleSeqAlign  MultipleSeqAlign::dropClones(){
	MultipleSeqAlign nmsa;
	string msaname="Noclone"+getName();
	nmsa.setName(msaname);
	map<string,Sequence> non_redundant_list; //Non-redundant list of top-id sequences found so far
	int oldsz=0;
	for(MSA::iterator st=beginSeq();st!=endSeq();st++){
		non_redundant_list[st->alignedSequence()]=*st;
		if(oldsz<non_redundant_list.size()){
			oldsz=non_redundant_list.size();
			nmsa.addSeq(*st);
		}
	}
	if(nmsa.empty()&&!QUIET)cout<<"#WARNING : NO sequence extracted"<<endl;
	return nmsa;
}

///For each sequence of MSA-B, choose among top 10 highest ID hits against MSA-A, the top one that hasn't been chosen yet.
MultipleSeqAlign  MultipleSeqAlign::xtractSequencesHighestId(MultipleSeqAlign* msab, int cullsize, vector<int>* positions){
	MultipleSeqAlign nmsa;
	string seqstr;
	nmsa.setName("Extracted sequences");
	bool cond,identical;
	double identicalPID=_Seqlist.begin()->id(*_Seqlist.begin()); //PID for exactly identical secuences
	double sid;
	double minId;
	double thrId=0.0;
	int oldsz=0;
	int seq_oldsz;
	int TOPLISTSIZE=100;
	if(VERBOSE)cout<<"#EXtracting top sequences for..."<<endl;
	map<string,Sequence> non_redundant_list; //Non-redundant list of top-id sequences found so far
	multimap<double,Sequence,greaterThan> topid; //Same as non_redundant_list but ordered by sequence ID starting from highest.
	for(MSA::iterator stb=msab->beginSeq();stb!=msab->endSeq();stb++){
		minId=0.0;
		seq_oldsz=0;
		map<string,Sequence> seq_non_redundant_list;
		multimap<double,Sequence,greaterThan> st_topid;
		if(VERBOSE)cout<<"# "<<stb->name()<<" : "<<endl;
		for(MSA::iterator st=_Seqlist.begin();st!=_Seqlist.end();st++){
			identical=false;
			for(map<string,Sequence>::iterator id=non_redundant_list.begin();id!=non_redundant_list.end();id++){
				if(st->id(id->second,positions)==identicalPID){; //Is contained in one of the already found sequences? If so
					identical=true;
					break;
				}
			}
			if(identical)continue;//... skip it as we want only NON-redundant
			seqstr=st->alignedSequence();
			sid=st->id(*stb,positions);
			if(VERBOSE)cout<<"#\t"<<st->name()<<"\t"<<sid<<endl;
			cond=(sid>minId);
			if(cond||seq_oldsz<TOPLISTSIZE){ //Cull top TOPLISTSIZE sequences from MSA-A
				if(seq_oldsz==TOPLISTSIZE){
					//seq_non_redundant_list.erase((seq_non_redundant_list.rbegin())->first); //remove last item
					seq_non_redundant_list.erase(st_topid.rbegin()->second.alignedSequence()); //remove last item
					seq_oldsz--;
					st_topid.erase((st_topid.rbegin())->first); //remove last item
				}
				seq_non_redundant_list[seqstr]=*st;
				if(seq_oldsz<seq_non_redundant_list.size()){ //Update lists if sequence hasn't been added before
					seq_oldsz=seq_non_redundant_list.size(); //update number of top-id sequences found so far. 
					st_topid.insert(pair<double,Sequence>(sid,*st)); //update top id list . We could use this map as well for updating oldsz
					minId=st_topid.rbegin()->first; //Update lowest top id value
					if(VERBOSE)cout<<"#Got : "<<st->name()<<" (seq_oldz="<<seq_oldsz<<") "<<"\t"<<sid<<"\t\t(minID="<<minId<<")"<<endl;
				}
			}
		}
		for(multimap<double,Sequence,greaterThan>::iterator it=st_topid.begin();it!=st_topid.end();it++){
			non_redundant_list[(it->second).alignedSequence()]=it->second; //Try to add sequence to non-redundant list
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
		seqstr=st->alignedSequence();
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
                                (*ita).printAlignment();
                                (*itb).printAlignment();
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
	bool sameseq=false;
        for(MSA::iterator ita=_Seqlist.begin();ita!=_Seqlist.end();ita++)
                for(MSA::iterator itb=msa.beginSeq();itb!=msa.endSeq();itb++){
                        string na=(*ita).name();
                        if(na.substr(0,1).compare(">")==0)na=na.substr(1,na.length()-1);
                        string nb=(*itb).name();
                        if(nb.substr(0,1).compare(">")==0)nb=nb.substr(1,nb.length()-1);
			sameseq=ita->alignedSequence().compare(itb->alignedSequence())==0?true:false;
                        cout<<na<<"\t"<<nb<<"\t"<<SeqId(*ita,*itb,positions)<<"\t"<<sameseq<<endl;
                }
}

void MultipleSeqAlign::printPairwiseIds(vector<int>* positions)
{
	bool sameseq=false;
        for(MSA::iterator ita=_Seqlist.begin();ita!=_Seqlist.end()-1;ita++)
                for(MSA::iterator itb=ita+1;itb!=_Seqlist.end();itb++){
                        string na=(*ita).name();
                        if(na.substr(0,1).compare(">")==0)na=na.substr(1,na.length()-1);
                        string nb=(*itb).name();
                        if(nb.substr(0,1).compare(">")==0)nb=nb.substr(1,nb.length()-1);
			sameseq=ita->alignedSequence().compare(itb->alignedSequence())==0?true:false;
                        cout<<na<<"\t"<<nb<<"\t"<<SeqId(*ita,*itb,positions)<<"\t"<<sameseq<<endl;
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
                (*ita).printAlignment();
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
