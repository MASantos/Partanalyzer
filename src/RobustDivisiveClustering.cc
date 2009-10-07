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

RobustDivisiveClustering::RobustDivisiveClustering(char* graphfile, bool below, int nsamples, pmetricv pmetric, double extensivity){
	_initialize(graphfile);
	runPruningClustering(below,nsamples, pmetric, extensivity);
}
RobustDivisiveClustering::RobustDivisiveClustering(MatrixOfValues* graph, bool below, int nsamples, pmetricv pm, double extensivity){
	_initialize(graph);
	runPruningClustering(below,nsamples, pm, extensivity);
}

void  RobustDivisiveClustering::_setMainParameters(bool prunebelow,int samples, pmetricv metric, double ext){
	below=prunebelow;
	nsamples=samples;
	pmetric=metric;
	extensivity=ext;
	///Step of pruning thresholds
	_delta=(_max-_min)/(1.0*nsamples);
	///Must be odd integer: 2x#neighbors+1
	_MOVINGWINDOW=2*numberOfNeighbors+1;
}

void RobustDivisiveClustering::_calculateVariation(vector<Partition >& windpart, int& refpart , double& gauge, bool reducewin){
	if(VERBOSE){
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
	if(VERBOSE){
		cout<<"#BestOptimal partition so far ("<<optimal_partitions.size()<<"):"<<endl;
		optimal_partitions.begin()->second.printPartition();
		cout<<"#-----------------------------------------"<<endl;
	}
	///FILO stack of partiions: constitute the averaging window of partitions
	if(windpart.size()==_MOVINGWINDOW||reducewin)windpart.erase(windpart.begin());
	//if(it>=_MOVINGWINDOW)windpart.erase(windpart.begin());
}

void RobustDivisiveClustering::runPruningClustering(bool prunebelow, int samples, pmetricv metric, double extensivity){
	_setMainParameters(prunebelow,samples,metric,extensivity);
	_summaryParameters();
	int refpart=0;
	///pruning treshold for refpart relative to current (last) partition
	double gauge=numberOfNeighbors*_delta;
	///Averaging window of partitions
	vector<Partition > windpart;
	///Perform nsamples steps with a pruning threshold > minimum edge value
	for(int it=0;it<=nsamples;it++){
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
			P = _graphptr->pruneEdgesBelow(_graph_pruning_thr, false).cluster();
		}else{
			///or pruning above threshold
			P = _graphptr->pruneEdgesAbove(_graph_pruning_thr, false).cluster();
		}
		///... and add it to the averaging window
		windpart.push_back( P );
		if(DEBUG){
			P.summary();
			P.printPartition();
		}
		///If we have at least numberOfNeighbors + 1, then we can do calculations
		if(it>=numberOfNeighbors){
			///Index of reference partition  (index starts by 0).
			refpart=(it+1>=_MOVINGWINDOW)?numberOfNeighbors:(it-numberOfNeighbors);
			if(VERBOSE)cout<<"#it="<<it<<" refpart="<<refpart<<" winsize="<<windpart.size()<<" #neigh="<<numberOfNeighbors<<" Movingwindow="<<_MOVINGWINDOW<<" nsamples="<<nsamples<<endl;
			///Calculate average distance between neighbors for partition refpart in window
			_calculateVariation(windpart, refpart, gauge);
		}
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
