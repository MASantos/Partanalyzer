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

#ifndef _CLASS_ROBUSTDIVISIVECLUSTERING_H
#define _CLASS_ROBUSTDIVISIVECLUSTERING_H 1

#include "PartitionStats.h"
#include "MatrixOfValues.h"
#include "Statistics.h"

#define RDC_DEFAULT_NUMBER_SAMPLES 100
#define RDC_DEFAULT_NUMBER_NEIGHBORS 2
#define RDC_DEFAULT_TOP_BEST_PARTITIONS 20
/**Performs a robust pruning clustering of the input graph
 * */
class RobustDivisiveClustering
{
	MatrixOfValues* _graphptr;
	bool _inStack;
	double _min;
	double _max;
	double _delta;
	double _variance;
	double _graph_pruning_thr;
	int  _MOVINGWINDOW;
	void _initialize(char* graphfile);
	void _initialize(MatrixOfValues* graph);
	void _setMainParameters(bool prunebelow,int samples, pmetricv metric, double ext, bool selfConsistently);
	void _summaryParameters();
	void _calculateVariation(vector<Partition >& windpart, int& refpart, double& gauge, bool reducewin=false);
	void runPruningClustering(bool below, int nsamples, pmetricv pm, double extensivity, bool selfConsistently);
public:
	RobustDivisiveClustering(char* graphfile, bool prunebelow, int nsamples=RDC_DEFAULT_NUMBER_SAMPLES, pmetricv pm=shannon, double extensivity=EXTENSIVITY_DEFAULT, bool selfConsistently=false);
	RobustDivisiveClustering(MatrixOfValues* graph, bool prunebelow, int nsamples=RDC_DEFAULT_NUMBER_SAMPLES, pmetricv pm=shannon, double extensivity=EXTENSIVITY_DEFAULT, bool selfConsistently=false);
	~RobustDivisiveClustering(); 
	//RobustDivisiveClustering(RobustDivisiveClustering& RDC);
	//RobustDivisiveClustering operator=(RobustDivisiveClustering& RDC);
	void reset();
	bool below;
	double extensivity;
	int nsamples;
	pmetricv pmetric;
	int numberOfNeighbors;
	int topBestPartitions;
	map<double, Sampling > pdistance_vs_pruningthr;
	map<double, Partition> optimal_partitions;
	int optimal_Nsamples(int Nbsamples, bool below);
	//map<double, double> optimal_pruningthr;
	void printPartition();
	void printDistanceVsPruningThreshold();
};

#endif //END _CLASS_ROBUSTDIVISIVECLUSTERING_H
