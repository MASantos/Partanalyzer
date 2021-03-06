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

#ifndef _CLASS_CCOP_H
#define _CLASS_CCOP_H 1

#include "partanalyzer_includes.h"
#include "partanalyzer_definitions.h"
#include "partanalyzer_basic_operations.h"

#include "Partition.h"
#include "MatrixOfValues.h"

///DEFAULT THERMODYNAMIC COEFF.
//extern double beta;
//extern double mu;

///Checks the consistency of a partition relative to a given graph (matrix of edge values).
/**Checks the consistency of a partition in relation to a given graph (matrix of values). It's the initial core of this whole project.
 * */
class ccop
{
	MatrixOfValues* _MX;
	Partition* _part;
	double _w_intra,_w_inter,_w_intra_thr,_f_intra_thr,_threshold;
	long int _nintra,_ninter;
	double _chi2;

	void _checkClusters(\
		long int& K, int& ofs,\
		long int& total_nodes, int& non_singletons,\
		long int& non_singleton_nodes, double& cw_intra,\
		double& cw_inter, double& cf_intra_thr, \
		double& cw_intra_thr, double& w_intra_thr, \
		double& cw_chi2, double& c_chi2, \
		long int& nintra, long int& ninter, \
		double& sw_inter, double& w_intra, \
		double& w_inter, double& tf_intra_thr, \
		double& pa, double& S, \
		double& E, double& Cv, \
		double& cE, double& N, \
		double& Z, double& Es, \
		double& Ei, double& expected, \
		double& c_Qval, double& Qval, \
		int& nstableclusters, long int& nstableclusters_size\
		);

	void _evaluateCluster(\
		smat::iterator cla , long int& K, \
		int& ofs, int& non_singletons,\
		long int& non_singleton_nodes,double& cw_intra,\
		double& cw_inter, double& cf_intra_thr, \
		double& cw_intra_thr, double& w_intra_thr, \
		double& cw_chi2, double& c_chi2, \
		long int& nintra, long int& ninter, \
		double& sw_inter, double& w_intra, \
		double& w_inter, double& tf_intra_thr, \
		double& pa, double& S, \
		double& E, double& Cv, \
		double& cE, double& N, \
		double& Z, double& Es, \
		double& Ei, double& expected, \
		double& c_Qval, double& Qval, \
		int& nstableclusters, long int& nstableclusters_size\
		);

	void _getIntraClusterValuesPerElement(\
		svect::iterator ita, smat::iterator cla, \
		long int& nintra, double& cw_intra, \
		double& cf_intra_thr, double& cw_intra_thr\
		);

	void _getInterClusterValuesPerElement(\
		svect::iterator ita, smat::iterator cla, \
		int& ofs, long int& ninter, \
		double& cw_inter, double& sw_inter, \
		double& pa, double& S, \
		double& E, double& Cv\
		);

public:
	ccop(MatrixOfValues* MX, Partition* part);
	ccop(MatrixOfValues* MX, Partition* part, double thr);
	void checkConsistency();
	
	void distribution();
	//double mx(string a, string b){ return _MX(a,b);}
	double threshold(){return _threshold;}
	//vector< double > observables;
	bool isSparseGraph;
};

#endif //END _CLASS_CCOP_H
