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
#include "partanalyzer_includes.h"
#include "partanalyzer_definitions.h"
#include "partanalyzer_basic_operations.h"

#ifndef _CLASS_PARTITION_H
#define _CLASS_PARTITION_H 1


/**Main class Partition implements a partition and its algebra
 * */
class Partition
{
	///Filename containing partition
	char* _partitionf;
	//PART Partition format offset
	int _items_offset;
	///Number of clusters
	long int _nclusters;
	///Number of singletons
	long int _nsingletons;
	///Number of clusters with only 2 elements
	long int _npairs;
	/**Number of non-trivial clusters (Defined here as those with 4 or more elements. 
	///This assumes that the majority of the pairs can be retrieved using any reasonable tree-building algorithm by
	///cutting just above the first level nodes. Maybe even some of the triples, by judicious selection of isolated branch leafs 
	///from next level nodes. However, is seems much less likely to recover all (many) clusters with 4+ elements unless the
	///algorithm is really fine tuned.
	*/
	long int _nnontrivial;
	///Number of elements (of underlying space) 
	long int _nitems;
	///Index of largest cluster in member smat clusters
	int _largest_cluster;
	///Iterator pointing to largest cluster in member smat clusters
	smat::iterator _it_largest_cluster;
	///Associated Adjacency matrix
	graph _Ad;
	///Update all member elements
	void _resetMembers();
	///Read partition from a file
	void _readClusters();
	///Phony function depending on file format; Currently, print read format if VERBOSE
	void partitionInputFormat(partFileFormat iformat);
	svect _it_found;
	///Input format for reading partitions: MCL (partFmtMCL), FREE (partFmtFREE) and the default one, partanalyzer's own format (partFmtPART)
	partFileFormat _piformat;
	///Specific for MCL format: the partition tab file maps element labels to integers
	char* _mcltabf;	
	///The MCL tab file
	map<int, string> _mcltab;
	///Unordered list of items
	vector<string > _vitems;
public:
	///Vector of clusters, each being a svect containing number of elements, name of cluster, and all its elements.
	smat clusters;
	///Set of elements (treated as strings; represents the underlying set of elements).
	sset sitems;
	///Sorted set of singleton elements.
	sset ssingletons; ///Let's make sure we do not get any duplicated singleton. This is just a cheap patch for trying to speed up isaPartition().
	///Default Partition constructor
	Partition();
	///Instantiate partition from file using format iformat (uses default offset for partanalyzer's own format partFmtPART)
	// Partition(char* file, partFileFormat iformat);
	///Instantiate partition from file using format iformat and specified offset (applies only to partanalyzer's own format partFmtPART)
	Partition(char* file, partFileFormat iformat=partFmtPART, int ofs=2);
	///Create a partition out of a set of clusters
	Partition(smat* clustersl, int ofs, bool dosort=true, char* partf=NULL , char* tabf=NULL);
	///Extract selected elements
	Partition xtractElements(svect* elements);
	///Create a consensus partition 
	void xtrConsPart(multimap<int,string,greaterThan> consPart, int ofs=1);
	///Get (and reset) the underlying set of elements sorted by alphanumerically
	sset getItems();
	///Get an ordered set with all elements contained in (svect) cl. {Doesn't affect Partition; should be external maybe}
	sset getItems(cluster* cl);
	///Check it is a sound Partition of the specified set of elements
	bool isaPartitionOf(svect& ssetOfElements); 
	///Check it is a sound Partition of the default set of elements
	bool isaPartition();
	///Print partition using the default partanalyzer format
	void printPartition(bool SequentialClusterNames=false, string ClusterPrefix="C");
	///Print partition in the specified format 
	void printPartition(partFileFormat iformat, bool SequentialClusterNames=false, string ClusterPrefix="C");
	///Get number of clusters (including singletons)
	long int n_clusters(){ return _nclusters;}
	///Get number of singletons
	long int n_singletons(){ return _nsingletons;}
	///Get number of pair clusters
	long int n_pairs(){ return _npairs;}
	///Get number of non-singletons clusters
	long int n_nonSingClusters(){ return _nclusters-_nsingletons;}
	///Get number of elements
	long int n_items(){ return _nitems;}
	///Get cluster offset (used only with partanalyzer's own format partFmtPART)
	int cluster_offset(){return _items_offset;}
	///Get index of largest cluster in (smat) Partition::smat
	int largest_Cluster(){ return _largest_cluster;}
	///Get iterator pointing to largest cluster in (smat) Partition::smat
	smat::iterator it_largest_Cluster(){ return _it_largest_cluster;}
	///Set partition input file format
	void setPartInputFormat(partFileFormat iformat){ _piformat=iformat;}
	///Retrieve partition input file format
	partFileFormat getPartInputFormat(partFileFormat iformat){ return _piformat;}
	///Load MCL tab file
	void mclTabFile(char* mcltabf);
	///Swap labels according to the provided tab file
	void swapLabels(char* mcltabf);
	///Entropy of a partition
	double H();
	///Cardinality of a partition
	long int card(){return n_clusters();}
	///Boltzman Entropy: log (#available states)
	double BK(){return log( (double)1.0*card() ) ;}
	///Tsallis Entropy: is equal to H iif q=1 ; q=non-extensivity coefficient
	double TS(double q=EXTENSIVITY_DEFAULT_TSALLIS);
	///Renyi Entropy: is equal to H iif q=1 ; q=non-extensivity coefficient
	double RS(double q=EXTENSIVITY_DEFAULT_RENYI);
	///Calculates the Jeffrey Qnorm (a la Tarantola) of a partition= exp(RS(q)) or expr(TS(q))
	double JQnorm(double q);
	///Calculates and prints vi distance againts given partition part2
	void vipp(Partition* part2);
	///Calculates and prints edit score distance againts given partition part2
	void edsc(Partition* part2);
	///Calculates the Tarantola distance against given partition part2
	void td(Partition* part2, double q);
	///Calculates the purity strict values againts the given partition p2
	double purityStrict(Partition* part2){ return purityScore(part2)[0];}
	///Calculates the purity lax values againts the given partition p2
	double purityLax(Partition* part2){ return purityScore(part2)[1];}
	///Calculates both purity score values againts the given partition p2
	vector<double> purityScore(Partition* part2);
	void missing(svect* items_found);
	void missing();
	///Given an element, get the size of the cluster that contains it.
	int getClusterSize(const string& item);
	///Given an element, get the name of the cluster that contains it.
	string getClusterName(string& item);
	///Given an index from clusters, get the name of that cluster
	string getClusterName(int& clidx);
	///Given the cluster, get its name
	string getClusterName(svect& cluster);
	///Given an element, get the index of the cluster that contains it.
	int getClusterIdx(string& item);
	///Two elements are equivalent if they belong to the same cluster
	bool areEquiv(string a, string b);
	///Return name of cluster if true; x if false ; NAN1 (NAN2) if first (second) not found.
	string areWithinSameCluster(string ita, string itb);
	///Get the filename.
	char* FileName(){ return _partitionf;}
	///cluster is a typedef specific for Partition
	svect getClusterOf(string item);
	///Build associated Adjacency matrix
	graph setAdjacencyMatrix();
	///Build associated Adjacency matrix sorting elements by clusters
	graph setAdjacencyMatrix_os();
	///Get AdjacencyMatrix
	graph Ad(){ if(_Ad.empty()) return setAdjacencyMatrix(); return _Ad;}
	///Get AdjacencyMatrix with indexes ordered as given by the partition
	graph Ad_os(){ if(_Ad.empty()) return setAdjacencyMatrix_os(); return _Ad;}
	//
	// RELATIONS AND OPERATIONS ON PARTITOINS 
	//
	///Project a Partition of N elements onto a subspace of M<=N elements (active view)
	void SubsProject(sset& itemset);
	///Calculate the partition intersection of the present one and part2 (passed by pointer)
	Partition intersection(Partition* part2);
	///Calculate the partition intersection of the present one and part2 (passed by reference)
	Partition intersection(Partition& part2){ Partition z=part2; return intersection(&z); }
	///Define Lattice preorder relation (as an interface)
	bool lessThan(Partition* part2){ return purityLax(part2)==1?true:false;}
	///Define Lattice preorder relation (as a binary operator)
	bool operator<=(Partition& part2);
	///Define lattice strict preorder relation (as a binary operator)
	bool operator<(Partition& part2){ return (n_clusters()>part2.n_clusters()&& *this<=part2)?true:false;}
	///Define equality operater
	bool operator==(Partition& part2);
	///Define intersetion of partitions ( as a binary operator ; arguments as refereces)
	Partition operator*(Partition& part2){ return intersection(&part2); }
	///Define intersetion of partitions ( as a binary operator ; arguments as pointers)
	Partition operator*(Partition* part2){ Partition& p=*part2; return operator*(p);}
	///(TO BE IMPLEMENTED) Define union of partitions (as a binary operator)
	Partition operator+(Partition& part2);
};
#endif //END  _CLASS_PARTITION_H

