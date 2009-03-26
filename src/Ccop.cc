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

#ifndef _CLASS_CCOP
#define _CLASS_CCOP 1

#include "Ccop.h"

/**Checks the consistency of a partition in relation to a given graph (matrix of values). It's the initial core of this while project.
 * */
ccop::ccop(MatrixOfValues* MX, Partition* part){
	_MX=MX;
	_part=part;
	_threshold=-1.0;
	_chi2=0.0;
	long int posNedges2=_part->n_items()*(_part->n_items()-1);
	isSparseGraph=(2*_MX->n_edges()< posNedges2)?true:false;
	//_runCheck();
}

ccop::ccop(MatrixOfValues* MX, Partition* part, double thr){
	_MX=MX;
	_part=part;
	_threshold=thr;
	_chi2=0.0;
	long int posNedges2=_part->n_items()*(_part->n_items()-1);
	isSparseGraph=(2*_MX->n_edges()< posNedges2)?true:false;
	//_runCheck();
}

void ccop::distribution(){
	int ofs,cln;
	ofs=_part->cluster_offset();
	cout<<"#BeginIntraClusterDistribution"<<endl;
	for(smat::iterator cla=_part->clusters.begin();cla!=_part->clusters.end();cla++)
		if(cla->size()>1+ofs) //Ignore singletons
			for(svect::iterator ita=cla->begin()+ofs;ita!=cla->end();ita++)
				for(svect::iterator itb=ita+1;itb!=cla->end();itb++)
					cout<<*ita<<"\t"<<*itb<<"\t"<<_MX->v(*ita,*itb)<<endl;
	cout<<"#EndIntraClusterDistribution"<<endl;
	cout<<"#BeginInterClusterDistribution"<<endl;
	for(smat::iterator cla=_part->clusters.begin();cla!=_part->clusters.end();cla++)
		for(smat::iterator clb=_part->clusters.begin();clb!=_part->clusters.end();clb++)
			if(clb!=cla)
				for(svect::iterator ita=cla->begin()+ofs;ita!=cla->end();ita++)
					for(svect::iterator itb=clb->begin()+ofs;itb!=clb->end();itb++)
						cout<<*ita<<"\t"<<*itb<<"\t"<<_MX->v(*ita,*itb)<<endl;
	cout<<"#EndInterClusterDistribution"<<endl;
}

//void ccop::_runCheck(){
void ccop::checkConsistency(){
	int ofs,cln;
	int nstableclusters=0;
	long int nstableclusters_size=0;
	_w_intra=_w_inter=_f_intra_thr=0.0;
	_nintra=_ninter=0;
	long int nintra,ninter;
	double w_intra,w_inter,sw_inter,w_intra_thr,cw_intra_thr,tf_intra_thr;
	double cw_intra,cw_inter,cf_intra_thr;
	double cw_chi2,c_chi2,w_chi2,c_Qval,Qval;
	Qval=c_Qval=w_chi2=0.0;
	ofs=_part->cluster_offset();
	long int K=_part->clusters.size(); //number of clusters
	long int non_singleton_nodes,total_nodes;
	int non_singletons;
	non_singletons=0;
	non_singleton_nodes=total_nodes=0;
	cln=0;
	double expected;
	w_intra=w_inter=sw_inter=w_intra_thr=tf_intra_thr=0.0;
	double Z=0.0; // PARTITION FUNCTION
	double Cv=0.0;
	double cE,E,S;
	cE=E=S=0.0;
	double N=0.0;
	double pa;
	double Es,TE,Ei;
	Es=TE=Ei=0.0;
	cout.setf(ios::fixed,ios::floatfield);
	if(!QUIET)cout<<"#Calculating consistency check..."<<endl;
	cout<<"#partition-size= "<<_part->clusters.size()<<" offset= "<<ofs<<" threshold= "<<_threshold<<endl;
	//cout<<"#\n#Clusters averages: w_intra/w_inter\tw_intra\tw_inter\tintra_thr\tnintra\tninter"<<endl;
	//cout<<"#\n#Clusters averages: chi2\tw_chi2\tw_intra\tw_inter\tintra_thr\tnintra\tninter"<<endl;
	//cout<<"#\n#Clusters averages: chi2\tw_chi2\tintra\tinter\tw_intra\tw_inter\tintra_thr\tnintra\tninter"<<endl;
	cout<<"#\n#Clusters values: \n#index\tClname\tsize\tw_intra\t\tw_inter\t\t<w_intra>\t<w_inter>  f_intra_thr\tnintra\tninter"<<endl;
	for(smat::iterator cla=_part->clusters.begin();cla!=_part->clusters.end();cla++){
		total_nodes+=cla->size()-ofs;
		if(cla->size()>1+ofs){ //Do not check consistency of singletons
			non_singletons++;
			non_singleton_nodes+=cla->size()-ofs;
			cw_intra=cw_inter=cf_intra_thr=cw_intra_thr=0.0;
			cw_chi2=c_chi2=0.0;
			nintra=ninter=0;
			for(svect::iterator ita=cla->begin()+ofs;ita!=cla->end();ita++){
				for(svect::iterator itb=ita+1;itb!=cla->end();itb++)
					if(_MX->existEdge(*ita,*itb)){
						//Number of intra-cluster edges
						nintra++;
						//Intra-cluster weight 
						cw_intra+=_MX->v(*ita,*itb);
						if(_MX->v(*ita,*itb)>_threshold) { 
							//Fraction of intra-cluster edges above threshold
							cf_intra_thr+=1.0; 
							//Intra-cluster weight above threshold
							cw_intra_thr+=_MX->v(*ita,*itb) ;
						}
						if(VERBOSE) cout<<"#"<<*ita<<" - "<<*itb<<" = "<<_MX->v(*ita,*itb)<<endl;
					}
				for(smat::iterator clb=_part->clusters.begin();clb!=_part->clusters.end();clb++)
					if(clb!=cla)
						for(svect::iterator itb=clb->begin()+ofs;itb!=clb->end();itb++)
							if(_MX->existEdge(*ita,*itb)){ 
								//Number of inter-cluster edges
								ninter++;
								//Inter-cluster weight 
								cw_inter+=_MX->v(*ita,*itb);
								if(clb->size()==1+ofs) {
									//Inter-cluster weight between a non-singleton and a singleton
									sw_inter+=_MX->v(*ita,*itb);
								//Thermodynamic magnitudes
									//Ensemble probability
									pa=exp(-beta*(-0.5*_MX->v(*ita,*itb)-mu));
									//Entropy
									S+=-pa*log(pa);
									//Internal energy
									E+=-0.5*_MX->v(*ita,*itb)*pa;
									//Heat capacity
									Cv+=0.25*_MX->v(*ita,*itb)*_MX->v(*ita,*itb)*pa;
									//Chi-square
									_chi2+=pow(1.0-_part->n_items()*pa,2);
								}
							}
			}
			cE=-cw_intra/K-0.5*cw_inter;
			//double pa=(1.0*cla->size()-ofs)/_part->n_items(); // probability of cluster
			pa=exp(-beta*(cE-mu*(1.0*cla->size()-ofs)));
			Z+=pa;
			S+=-pa*log(pa);
			/// Thermodynamic averages
			//N+=(1.0*cla->size()-ofs)*pa;
			N+=_part->n_items()*pa;
			w_intra+=cw_intra*pa;
			E+=cE*pa;
			Cv+=cE*cE*pa;
			w_inter+=0.5*cw_inter*pa;
			w_intra_thr+=cw_intra_thr*pa;
			tf_intra_thr+=(cf_intra_thr/nintra*100.0)*pa;
			//Surface tension
			Es+=0.5*cw_intra;
			Ei+=cw_inter;
			///Expected (random) in/out-going weight
			//expected=pow(cw_intra+cw_inter,2); 
			//cw_chi2=pow(cw_intra-expected,2)/expected;
			expected=1.0*cla->size()-ofs;
			c_chi2=pow(_part->n_items()*pa-expected,2)/expected;
			///Expected (random) in/out-going edges
			expected=(nintra+ninter)==0?0:pow((double)(nintra+ninter),2);
			//c_chi2=pow(nintra-expected,2)/expected;
			c_Qval+=expected;
			/// Newman & et al. modularity value
			Qval+=(double)nintra;
			/// Cluster averages
			if(nintra>0)cw_intra/=(double)nintra;
			if(ninter>0)cw_inter/=(double)ninter;
			if(nintra>0)cf_intra_thr/=(double)nintra/100.0;
			//c_chi2 becomes inf too often...
			//cout<<non_singletons<<"\t"<<(*cla)[1]<<"\t"<<cla->size()-ofs<<"\t"<<cw_intra*nintra<<"\t"<<cw_inter*ninter<<"\t"<<cw_intra<<"\t"<<cw_inter<<"\t"<<(int)cf_intra_thr<<"\t"<<nintra<<"\t"<<ninter<<"\t"<<c_chi2<<endl;
			cout<<non_singletons<<"\t"<<(*cla)[1]<<"\t"<<cla->size()-ofs<<"\t"<<cw_intra*nintra<<"\t"<<cw_inter*ninter<<"\t"<<cw_intra<<"\t"<<cw_inter<<"\t"<<(int)cf_intra_thr<<"\t"<<nintra<<"\t"<<ninter<<endl;
			//w_chi2+=cw_chi2;
			_chi2+=c_chi2;
			// Overall Cluster averages
			_w_intra+=cw_intra;
			_w_inter+=cw_inter;
			_f_intra_thr+=cf_intra_thr;
			//Don't know what was this meant for...
			if(cw_inter==0){
				nstableclusters++;
				nstableclusters_size+=cla->size()-ofs;
			}
		}
	}
	TE=-(Es+Ei+sw_inter);
	Es=-(Es-Ei-sw_inter);
	double zeta=Es/TE;
	Z+=_part->n_singletons();
	N+=_part->n_items()*_part->n_singletons();
	S=S/Z+log(Z);
	E/=Z;
	Cv=Cv/Z-(E*E);
	w_inter+=0.5*sw_inter;
	if(!QUIET)cout<<"#TEST: nsing="<<_part->n_singletons()<<" nitems="<<_part->n_items()<<" sw_inter="<<sw_inter<<endl;
	///expected is a dummy variable, so we can use it here too and avoid declaring more variables.
	//expected=0.5*non_singleton_nodes*(non_singleton_nodes-1.0);
	expected=0.5*total_nodes*(total_nodes-1.0);
	c_Qval/=pow(expected,2);
	Qval/=expected;
	Qval-=c_Qval;
	cout<<"#Thermodynamic averages: Beta= "<<beta<<" mu= "<<mu<<" F="<<beta*E-S-mu*N/Z<<" logZ= "<<log(Z)<<endl;
	cout<<"#w_intra\tw_inter\tw_intra_thr\tg=w_intra_thr/w_intra\tf_intra_thr\tN\tCv\tE"<<endl;
	cout<<w_intra/Z<<"\t"<<w_inter/Z<<"\t"<<w_intra_thr/Z<<"\t"<<w_intra_thr/w_intra<<"\t"<<tf_intra_thr/Z<<"\t"<<N/Z<<"\t"<<Cv<<"\t"<<E<<endl;
	cout<<"#Partition averages: "<<endl;
	cout<<"#TE\tEs\tZeta\t\tQval\tw_intra\t\tw_inter\tf_intra_thr\tnon_singletons\t#non-singleton-nodes\tfraction-non-singleton-nodes\t#total-nodes"<<endl;
	cout<<TE<<"\t"<<Es<<"\t"<<zeta<<"\t"<<Qval<<"\t"<<_w_intra/non_singletons<<"\t"<<_w_inter/non_singletons<<"\t"<<_f_intra_thr/non_singletons<<"\t"<<non_singletons<<"\t"<<non_singleton_nodes<<"\t"<<non_singleton_nodes*1.0/total_nodes<<"\t"<<total_nodes<<endl;
	//cout<<"#stable_clusters/size= "<<nstableclusters<<" "<<nstableclusters_size<<"\t"<<_chi2<<endl;
	cout<<-Ei<<"\t"<<-sw_inter<<endl;
}

#endif //END _CLASS_CCOP
