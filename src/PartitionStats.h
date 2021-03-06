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


#include <algorithm>

#ifndef _CLASS_PARTITIONSTATS_H
#define _CLASS_PARTITIONSTATS_H 1

#include "partanalyzer_includes.h"
#include "partanalyzer_definitions.h"
#include "partanalyzer_basic_operations.h"

#include "Partition.h"
#include "sNeighborhood.h"
#include "BellNumber.h"

#include "Statistics.h"
#include "MatrixOfValues.h"

// Prints system date; is defined elsewhere : partanalyzer_help.*
extern void systemDate();

///Type definition of a vector of partitions. Useful with trying to process multiple partitions at once, like in class PartitionStats
typedef vector<Partition > ppvect;

///Implements analysis and operations between multiple partitions. 
/**Implements operations between multiple partitions. 
Calculates different kind of pair-wise potentials, distances, inequality-residuals, the cover, ...
*/
class PartitionStats
{
	///List of file that contain our partitions
	vector<Charr > _fnamel;
	///Main member: Vector of (instantiated) Partitions
	vector<Partition > _partitionl;
	///2nd List of file that contain paritions: used for comparing two sets of partitions
	vector<Charr > _fnamel2;
	///2nd Vector of (instantiated) Partitions: used for comparing two sets of partitions
	vector<Partition > _partitionl2;
	///Cover of underlying set: Maps a sNeighborhood to each element
	map<string,sNeighborhood> _cover;
	///The consensus partition in raw format (offset=1) order by cluster size from largest to smallest.
	multimap<int,string,greaterThan> _consPart;
	//Partition _consensusPartition;
	///Hasse diagram/hierarchy of partitions
	map<long int, ppvect > _hasseNodes;
	///Consensus Adjacency matrix (graph): Assumes all partitions defined on the same underlying set
	graph _Ad;
	///External matrix of Weights. Allows to define weighted potentials/probabilities.
	MatrixOfValues* _W;
	///Number of partitions. Equals _partitionl.size()
	int _npart;
	///General cluster offset
	int _clsnofs;
	///Tsallis/Renyi Extensivity degree
	double _extensivity_degree;
	//
	bool SETUPCONSENSUSP;
	//
	bool _DIST_SUBSPROJECT;
	///Resets to defaults, erase all content
	void _resetPartitionMembers();
	//
	long int _pmetric(Partition& p1, Partition& p2, const int f);
	double _pmetric(Partition& p1, Partition& p2, const double f);
	double _pmetric(Partition& p1, Partition& p2, pmetricv metric);
	double _f(Partition& p, pmetricv metric);
	//
	void _fullConstructor(vector<Charr > fnames, partFileFormat iformat=partFmtPART, double extensivity=EXTENSIVITY_DEFAULT, int ofs=CLUSTEROFFSET_DEFAULT, int clstat_normalization_ofs=0, char* mcltabfile=NULL);
	//
	void _fillPartitionList(vector<Charr > fnames, vector<Partition >& partlist, char* mcltabfile=NULL, partFileFormat iformat=partFmtPART, int ofs=CLUSTEROFFSET_DEFAULT);
	///Print header info lines before printing distance results
	void distancesPrintHeadComment(pmetricv metric, flagheader hd, bool usingREF=false);
	//
	void distancesPrintHeadComment(pmeasure measure, pmetricv metric, flagheader hd, bool usingREF=false);
public:
	BellNumber BellN;
	///Default constructor simply resets all members to "zero"
	PartitionStats();
	///Instantiates a PartitionStats out of a list of filenames, a common file input format, a float parameter, a common cluster-offset value and a common normalization factor gauge (default=0).
	PartitionStats(vector<Charr > fnames, partFileFormat iformat=partFmtPART, double extensivity=EXTENSIVITY_DEFAULT, int ofs=CLUSTEROFFSET_DEFAULT, int clstat_normalization_ofs=0, char* mcltabfile=NULL);
	///Instantiates a PartitionStats out of a list of filenames, a common file tabfile, and a float parameter used for non-standard entropies
	PartitionStats(vector<Charr > fnames, char* mcltabfile=NULL, double extensivity=EXTENSIVITY_DEFAULT);
	///Instantiates from a vector of partitions
	PartitionStats(vector<Partition >& vpart, double extensivity=EXTENSIVITY_DEFAULT);
	///Declares a second list of input partition files, or the initial one, if it hasn't been declared before
	void defineListOfPartitions(vector<Charr > fnames, char* mcltabfile=NULL, partFileFormat iformat=partFmtPART, int ofs=CLUSTEROFFSET_DEFAULT);
	///Erases second list of input partition files
	bool clear2ndListOfPartitions();
	///Checks whether a second list of input partitions file has been defined
	bool has2ndListOfPartitions(){ return (_partitionl2.size()>0);}
	///Checks if each of the provided partitions is a sound partition, i.e., if all clusters are pair-wise disjoint.
	int arePartitions();
	///Builds the Meet Partition (overall intersection) of the underlying set of elements out of the list of partitions.
	Partition getMeet();
	///Builds the Join Partition (overall union) of the underlying set of elements out of the list of partitions.
	Partition getJoin();
	///Builds the cover Partition (overall union) of the underlying set of elements out of the list of partitions.
	Partition getCoverPartition(){ return getJoin();}
	///Builds the cover _cover of the underlying set of elements out of the list of partitions.
	void getCover();
	///Check if _cover has been build.
	bool hasaCover(){ if(_cover.size()>1)return true; return false ;}
	///Prints out _cover, eventually printing out also the consensus if SETUPCONSENSUSP=true within the same loop. Notice this doesn't instantiates a consensus partition, but just prints it.
	void printCover(bool SETUPCONSENSUSP);
	///Used by printCover to print the consensus partition 
	void printConsensusPart();
	//void getConsensusPartition();
	///Explicitly build and instantiate a new partition which is the consensus one of the given list of partitions. 
	Partition getConsensusPartition();
	///Build consensus adjacency matrix
	void setConsensus_Ad();
	///Returns consensus adjacency matrix
	graph Ad(){ if(_Ad.empty()) setConsensus_Ad(); return _Ad;}
	///Print grey-scale image (pgm) of consensus adjacency matrix
	void pgm_Ad(){cout<<"#BeginConsensusAdjacencyMatrixPGMImage\n"<<_Ad<<"#EndConsensusAdjacencyMatrixPGMImage"<<endl;}
	///Print consensus adjacency matrix in raw format: list of rows for each edge : strA strB double
	void get_Ad();
	///Print fuzzy partition associated to the consensus adjacency matrix
	void get_FuzzyConsensusPartition();
	///Calculate cosine distance between associated Adjacency matrices
	void getAdCos();
	///Calculate all pair-wise purity scores
	void getPurity();
	///Calculate purity scores of all againts the first partition (reference)
	void getPurityRef();
	///Calculate purity scores of this partition (target) againts all other partition
	void getPurityTarget();
	///Compare cluster by cluster against the first partition's clusters (reference) and print overlap (fraction elements present in referece cluster)
	void getSplitsRef(splitmethod similarity);
	///get cardinality of partition
	long int card(Partition p){ return p.card();}
	///Get Shannon entropy of partition
	double H(Partition p){ return p.H();}
	///Get Boltzman entropy of partition
	double BK(Partition p){ return p.BK();}
	///Get Tsallis entropy of partition
	double TS(Partition p, double q){ return p.TS(q);}
	double TS(Partition p){ return p.TS();}
	///Get Renyi entropy of partition
	double RS(Partition p, double q){ return p.RS(q);}
	double RS(Partition p){ return p.RS();}
	///Get Jeffreys Qnorm (a la Tarantola) of partition
	double JQnorm(Partition p, double q){ return p.JQnorm(q);}
	///Print the value of the (information theoretic) potential associated with each partition
	void iPotential(pmetricv pm);
	///Calculates edit score distance between part1 and part2 using pmetric.
	long int ES(Partition& p1, Partition& p2){ return _pmetric(p1,p2,0);}
	///Calculates VI (Shannon) distance between part1 and part2 using pmetric. Should give the same as using partition explicitly built-in vipp function. Already checked?
	//double VI(Partition& p1, Partition& p2){return _pmetric(p1,p2,0.0);}
	double VI(Partition& p1, Partition& p2){return _pmetric(p1,p2,shannon);}
	///Calculates Boltzman distance between part1 and part2 using pmetric.
	double BK(Partition& p1, Partition& p2){ return _pmetric(p1,p2, boltzmann);}
	///Calculates Tsallis distance between part1 and part2 using pmetric.
	double TS(Partition& p1, Partition& p2){ return _pmetric(p1,p2, tsallis);}
	///Calculates Renyi distance between part1 and part2 using pmetric.
	double RS(Partition& p1, Partition& p2){ return _pmetric(p1,p2, renyi);}
	///Calculates Tarantola distance between part1 and part2 using Jeffrey's Qnorm.
	double TD(Partition& p1, Partition& p2){ return _pmetric(p1,p2, jeffreyQnorm);}
	///Calculates all distances againts the specified reference partition. Argument pmetricv specifies which metric to use (VI, Edit score,...)
	void distances(Partition& p, pmetricv pm=shannon);
	///Calculates all distances againts the reference partition. This is the first partition read. Argument pmetricv specifies which metric to use (VI, Edit score,...)
	void distancesRef(pmetricv pm=shannon);
	///Calculates all pair-wise distances. Argument pmetricv specifies which metric to use (VI, Edit score,...)
	void distances(pmetricv pm=shannon);
	/**Meadian, avg, std, std.err, var, min , max and sample sizeof the pair-wise distances among its partitions.
	If a non negative index is provided, distances are calculated with respect to partition refindex
	*/
	Sampling distancesDistribution(pmetricv pm=shannon,int refindex=-1);

	///Aproximate calculation of pair-wise distance between paritions with different number of elements
	/// It stripps of all elements that aren't share and calculates the distance using simply the rest.
	void distances_Subsprojection(pmetricv pm=shannon); 
	///Equivalent one for distances againts a common reference partition.
	void distancesRef_Subsprojection(pmetricv pm=shannon);
	///Calculates different symmetric and non-symmetric pair-wise measures
	void pmeasures(pmeasure measure, pmetricv pm=shannon); 
	///Calculates different symmetric and non-symmetric pair-wise measures against a common reference partition
	void pmeasuresRef(pmeasure measure, pmetricv pm=shannon); 
	///Prints the Hasse diagram corresponding to the given list of partitions.
	void printHasseDiagram();
	///Print 
	void printHasseNodes();
};

#endif //END  _CLASS_PARTITIONSTATS_H
