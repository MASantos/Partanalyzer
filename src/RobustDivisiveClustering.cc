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

#ifndef _CLASS_ROBUSTDIVISIVECLUSTERING
#define _CLASS_ROBUSTDIVISIVECLUSTERING 1

#include "RobustDivisiveClustering.h"

void RobustDivisiveClustering::_summaryParameters(){
	if(!QUIET){
		string prune=below?"Below":"Above";
		string metric;
		switch(pmetric){
			case shannon:
				metric="shannon"; break;
			case cardinality:
				metric="cardinality"; break;
			case renyi:
				metric="renyi"; break;
			case tsallis:
				metric="tsallis"; break;
			case jeffreyQnorm:
				metric="jeffreyQnorm"; break;
				
		}
		cout<<"#BeginRobustDivisiveClustering"<<endl;
		cout<<"#Pruning= "<<prune<<" Nsamples= "<<nsamples<<" metric= "<<metric<<" extensivity= "<<extensivity<<" neighbors= "<<numberOfNeighbors<<endl;
		cout<<"#min= "<<_min<<" max= "<<_max<<" Delta= "<<_delta<<endl;
		cout<<"#RecommendedMaxNsamples= "<<_variance<<" (variance)"<<endl;
	}
}

void RobustDivisiveClustering::reset(){
	_inStack=false;
	double extensivity=EXTENSIVITY_DEFAULT;
	int nsamples=RDC_DEFAULT_NUMBER_SAMPLES;
	pmetric=shannon;
	//pmetric=shannon;
	numberOfNeighbors=RDC_DEFAULT_NUMBER_NEIGHBORS;
	topBestPartitions=RDC_DEFAULT_TOP_BEST_PARTITIONS;
	pdistance_vs_pruningthr.clear();
	optimal_partitions.clear();
	///Get range of pruning threshold values
	Sampling graphEdgeDistribution=_graphptr->edgeDistributionStats();
	_min=graphEdgeDistribution.minimum_value();
	_max=graphEdgeDistribution.maximum_value();
	_variance=graphEdgeDistribution.variance();
}

void RobustDivisiveClustering::_initialize(char* graphfile){
	_graphptr= new MatrixOfValues (graphfile);
	reset();
	_inStack=true;
}

void RobustDivisiveClustering::_initialize(MatrixOfValues* graph){
	_graphptr=graph;
	reset();
}

RobustDivisiveClustering::~RobustDivisiveClustering(){
	reset();
	if(_inStack) delete _graphptr;
}

RobustDivisiveClustering::RobustDivisiveClustering(char* graphfile, bool below, int nsamples, pmetricv pmetric, double extensivity, bool selfConsistenly){
	_initialize(graphfile);
	runPruningClustering(below,nsamples, pmetric, extensivity, selfConsistenly);
}
RobustDivisiveClustering::RobustDivisiveClustering(MatrixOfValues* graph, bool below, int nsamples, pmetricv pm, double extensivity, bool selfConsistenly){
	_initialize(graph);
	runPruningClustering(below,nsamples, pm, extensivity, selfConsistenly);
}

int RobustDivisiveClustering::optimal_Nsamples(int samples, bool below){
	long int nprunedges;
	long int onprunedges;
	int lastN=1;
	int lastgoodN=1;
	bool lastoneOK=true;
	int DeltaN=samples;
	bool thisoneOK;
	double delta;
	set<int> badNs;
	if(!QUIET){
		cout<<"#Finding optimal nsamples: No="<<samples<<" _min="<<_min<<" _max="<<_max<<endl;
		systemDate();
	}
	while(true){
		thisoneOK=true;
		onprunedges=0;
		delta=(_max-_min)/(1.0*samples);
		int it;
		for(it=0;it<=samples;it++){
			nprunedges=0;
			_graph_pruning_thr=_min+it*delta;
			if(below){
				_graphptr->pruneEdgesBelow(_graph_pruning_thr,nprunedges, false).cluster();
			}else{
				_graphptr->pruneEdgesAbove(_graph_pruning_thr,nprunedges, false).cluster();
			}
			if(VERBOSE&&it%10==0)cout<<"it="<<it<<" ptr="<<_graph_pruning_thr<<" delta="<<delta<<" prunedges="<<nprunedges<<" Netpruned="<<nprunedges-onprunedges<<endl;
			if(it>=numberOfNeighbors && nprunedges-onprunedges==0){
				thisoneOK=false;
				if(!QUIET)cout<<"#Not good: it="<<it<<": "<<_graph_pruning_thr<<"/"<<nprunedges<<"("<<_graph_pruning_thr-delta<<"/"<<onprunedges<<")"<<endl;
				badNs.insert(samples);
				break;
			}
			onprunedges=nprunedges;
		}
		int newn;
		if((thisoneOK&&!lastoneOK)||(!thisoneOK&&lastoneOK)){
			if(thisoneOK){
				lastgoodN=samples;
				lastoneOK=true;
			}else{
				lastoneOK=false;
			}
			if(DeltaN==1||DeltaN==-1) {
				systemDate();
				return lastgoodN;
			}
			newn=(samples+lastN)/2;
		}
		else{
			if(thisoneOK){
				lastgoodN=samples;
				newn=samples+DeltaN;
			}else{
				if(lastgoodN==1){
					newn=samples-4*DeltaN/5;
				}else{
					newn=(samples+lastgoodN)/2;
				}
			}
		}
		DeltaN=newn-samples;
		lastN=samples;
		samples=newn;
		if(!QUIET)cout<<"#Trying now: N="<<samples<<" DeltaN="<<DeltaN<<endl;
	}
	cout<<"ERROR: RobustDivisiveClustering::optimal_Nsamples : Couldn't find optimal Nsamples"<<endl;
	exit(1);
}

void  RobustDivisiveClustering::_setMainParameters(bool prunebelow,int samples, pmetricv metric, double ext, bool selfConsistently){
	below=prunebelow;
	if(selfConsistently){
		nsamples=optimal_Nsamples(samples,prunebelow);
		if(!QUIET)cout<<"#OptimalNsamples= "<<nsamples<<endl;
	}else{
		nsamples=samples;
	}
	pmetric=metric;
	extensivity=ext;
	///Step of pruning thresholds
	_delta=(_max-_min)/(1.0*nsamples);
	///Must be odd integer: 2x#neighbors+1
	_MOVINGWINDOW=2*numberOfNeighbors+1;
}

void RobustDivisiveClustering::_calculateVariation(vector<Partition >& windpart, int& refpart , double& gauge, bool reducewin){
	if(DEBUG){
		cout<<"#Averaging over following "<<windpart.size()<<" partitions:"<<endl;
		for(vector<Partition >::iterator p=windpart.begin();p!=windpart.end();p++){
			p->summary();
			p->printPartition();
			//windpart.rbegin()->printPartition();
		}
	}
	///Instantiate the set of partitions within the averaging window
	PartitionStats pstat(windpart, extensivity);
	/**Ex., for numberOfNeighbors=2, calculate distances between the 3rd partition within the window
	and its neighbors and next-to-nearest neighbors
	*/
	Sampling pdistr=pstat.distancesDistribution(pmetric, refpart);
	if(DEBUG)cout<<"#DistanceDistribution: "<<_graph_pruning_thr-gauge<<"\t"<<pdistr<<endl;
	if(!QUIET){
		cout<<"#"<<_graph_pruning_thr-gauge<<"\t"<<pdistr.mean()<<"\t";
		windpart[refpart].summary();
	}
	///Store these sampling results vs pruning threshold
	pdistance_vs_pruningthr.insert( pair<double, Sampling > (_graph_pruning_thr-gauge,pdistr) );
	///Update list of optimal partitions and sort it from smallest mean distance to largest
	optimal_partitions.insert( pair<double, Partition> (pdistr.mean(), windpart[refpart]) );
	///Update list of optimal pruning threshold values
	//to be done 
	///Keep only the topBestPartitions best partitions 
	if(optimal_partitions.size()>topBestPartitions){
		if(VERBOSE)cout<<"#saved more than topBestPartitions ("<<topBestPartitions<<"): removing worst"<<endl;
		map<double, Partition>::iterator last=optimal_partitions.end();
		last--;
		///... by erasing the worst one.
		optimal_partitions.erase( last );
	}
	if(DEBUG){
		cout<<"#BestOptimal partition so far ("<<optimal_partitions.size()<<"):"<<endl;
		optimal_partitions.begin()->second.printPartition();
		cout<<"#-----------------------------------------"<<endl;
	}
	///FILO stack of partiions: constitute the averaging window of partitions
	if(windpart.size()==_MOVINGWINDOW||reducewin)windpart.erase(windpart.begin());
	//if(it>=_MOVINGWINDOW)windpart.erase(windpart.begin());
}

void RobustDivisiveClustering::runPruningClustering(bool prunebelow, int samples, pmetricv metric, double extensivity, bool selfConsistently){
	_setMainParameters(prunebelow,samples,metric,extensivity, selfConsistently);
	_summaryParameters();
	long int nprunedges=0;
	long int onprunedges=0;
	int refpart=0;
	///pruning treshold for refpart relative to current (last) partition
	double gauge=numberOfNeighbors*_delta;
	///Averaging window of partitions
	vector<Partition > windpart;
	///Perform nsamples steps with a pruning threshold > minimum edge value
	for(int it=0;it<=nsamples;it++){
		nprunedges=0;
		///
		///New step: update pruning threshold
		///Notice: both it and _graph_pruning_thr always point to the right-end of the sampling window
		_graph_pruning_thr=_min+it*_delta;
		//if(!QUIET)cout<<"#it="<<it<<"\tPruningThreshold= "<<_graph_pruning_thr<<endl;
		Partition P ;
		/**Prune graph accordingly and obtain partition by either
		*/
		if(below){
			///pruning below threshold
			P = _graphptr->pruneEdgesBelow(_graph_pruning_thr,nprunedges, false).cluster();
		}else{
			///or pruning above threshold
			P = _graphptr->pruneEdgesAbove(_graph_pruning_thr,nprunedges, false).cluster();
		}
		///... and add it to the averaging window
		windpart.push_back( P );
		if(DEBUG){
			P.summary();
			P.printPartition();
		}
		///If we have at least numberOfNeighbors + 1, then we can do calculations
		if(it>=numberOfNeighbors){
			//cout<<"pthr= "<<_graph_pruning_thr<<" REE-pruned= "<<nprunedges<<" ("<<onprunedges<<") "<<endl;
			if(nprunedges-onprunedges==0&&!QUIET)cout<<"#WARNING: NO additional edges deleted "<<nprunedges<<" pthr= "<<_graph_pruning_thr<<endl;
			///Index of reference partition  (index starts by 0).
			refpart=(it+1>=_MOVINGWINDOW)?numberOfNeighbors:(it-numberOfNeighbors);
			if(VERBOSE)cout<<"#it="<<it<<" pthr="<<_graph_pruning_thr<<" nprunedges="<<nprunedges<<" refpart="<<refpart<<" winsize="<<windpart.size()<<" #neigh="<<numberOfNeighbors<<" Movingwindow="<<_MOVINGWINDOW<<" nsamples="<<nsamples<<endl;
			///Calculate average distance between neighbors for partition refpart in window
			_calculateVariation(windpart, refpart, gauge);
		}
			onprunedges=nprunedges;
	}
	///Calculate variation for the last numberOfNeighbors partitions (They'll have less than 2*numberOfNeighbors 
	for(int it=0;it<numberOfNeighbors;it++){
		gauge=(numberOfNeighbors-1-it)*_delta;
		_calculateVariation(windpart, refpart, gauge,true);
	}
	if(!QUIET)cout<<"#EndRobustDivisiveClustering"<<endl;
}

void RobustDivisiveClustering::printPartition(){
	if(optimal_partitions.size()==0){
		cout<<"ERROR: RobustDivisiveClustering::printPartition() : optimal_partitions.size()==0 "<<endl;
		exit(1);
	}
	if(optimal_partitions.size()<topBestPartitions){
		if(!QUIET)cout<<"#WARNING: Number of optimal partitions "<<optimal_partitions.size()<<" < #topBestPartitions(="<<topBestPartitions<<")"<<endl;
	}
	//if(!QUIET)cout<<"#BeginRobustDivisiveClusteringPartition"<<endl;
	cout<<"#BeginRobustDivisiveClusteringPartition"<<endl;
	//cout<<"#pruningThreshold= "<<optimal_pruningthr.begin()->second;
	optimal_partitions.begin()->second.summary();
	optimal_partitions.begin()->second.printPartition();
	//if(!QUIET)cout<<"#EndRobustDivisiveClusteringPartition"<<endl;
	cout<<"#EndRobustDivisiveClusteringPartition"<<endl;
}

void  RobustDivisiveClustering::printDistanceVsPruningThreshold(){
	//if(!QUIET)cout<<"#BeginRobustDivisiveClusteringDistanceVsPruningThreshold"<<endl;
	//if(!QUIET)cout<<"#PruningThreshold\tMedian\tMean\tStd\tMerr\tVar\tMin\tMax\tTotNbNeighbors"<<endl;
	cout<<"#BeginRobustDivisiveClusteringDistanceVsPruningThreshold"<<endl;
	cout<<"#PruningThreshold\tMedian\tMean\tStd\tMerr\tVar\tMin\tMax\tTotNbNeighbors"<<endl;
	for(map<double, Sampling >::iterator it=pdistance_vs_pruningthr.begin(); it!=pdistance_vs_pruningthr.end(); it++){
		cout<<it->first<<"\t"<<it->second<<endl;
	}
	//if(!QUIET)cout<<"#EndRobustDivisiveClusteringDistanceVsPruningThreshold"<<endl;
	cout<<"#EndRobustDivisiveClusteringDistanceVsPruningThreshold"<<endl;
}
#endif //END _CLASS_ROBUSTDIVISIVECLUSTERING
