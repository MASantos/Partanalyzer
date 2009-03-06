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

#ifndef _CLASS_MATRIXOFVALUES_H
#define _CLASS_MATRIXOFVALUES_H 1

#include "Partition.h"


/**Allows dealing with graphs. A matrix of values here, in each line, specifies an edge value for a given pair of elements .
 * It also allows to merge two such graphs: The output is a pair of values -each from each of provided graphs- for each pair
 * of elements found. This operation has an alternative form where each of these 4-tuples get a label stating the cluster (subfamily)
 * name, if the edge is an intra-cluster edge; x if it's and inter-cluster edge. This requires to provide a partition. This is like
 * coloring the graph. Therefore, the command line option -color. The color can also be NAN1(NAN2) if the first(second) element
 * does not exist in the provided partition.
 * */
class MatrixOfValues
{
	char* _mxofvf;
	ifstream _is;
	svect _items;
	smap	_pairs;
	graph  _graph;
	row _mx;
	long int _nedges;
	long int _nitems;
	edge _Tweight;
	int _getIndexOfItem(string str);
public:
	MatrixOfValues();
	MatrixOfValues(char* file);
	///Get File name
	char* FileName(){return _mxofvf ;}
	///Read matrix of values from file
	void readMxValues();
	///Print matrix of values to standard output
	void printMatrix();
	///Get value of k-th edge
	double v(int k) { return _mx[k]; }
	///Get value of edge between i-th and j-th elements
	double v(int i, int j); 
	///Get value of edge spanned by nodes a and b
	double v(string a , string b, REDMxVal useRED=useOrgRED);
	//double v(string a , string b);
	///Merge two matrix of values
	void merge(MatrixOfValues* matrix2);
	///Merge two matrix of values and add cluster information for each edge
	void merge(MatrixOfValues* mx2, Partition* pt);
	///Cull edges specified in list matrix2
	void cull(MatrixOfValues* matrix2);
	///Cull edges spanned by nodes a and b
	strpair cullEdge(string a, string b);
	///Analyze cluster
	void clusteranalysis();
	///Get the number of edges present in the graph
	long int n_edges(){ return _nedges;}
	///Get the number of items spanning the actual graph
	long int n_items(){ return _nitems;}
	///Checks wether the edge was defined or not
	bool existEdge(string a, string b);
	///Gets total sum of edge weights
	edge W(){ return _Tweight;}
};

#endif //END _CLASS_MATRIXOFVALUES_H
