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

/** Partanalyzer
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2011.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )

Compile Options:

g++ -o partanalyzer partanalyzer.cc
*/

#ifndef _PARTANALYZER_DEFINITIONS_HEADER
#define _PARTANALYZER_DEFINITIONS_HEADER 1

#include "partanalyzer_includes.h"

///Current version
extern const char* VERSION;
///If true program outputs additional information regarding its calculations. This can be already quite a lot.
extern bool 	VERBOSE;
///If true program outputs even more detailed information
extern bool 	DEBUG;
///If true, program does not print any comment lines (those starting with #). No warnings will be printed either.
extern bool    QUIET;
///If true, we want to print some info that usually appears as a comment without the leading # comment sign
extern bool INFO;
///If true, we don't need to check whether any two clusters are mutually disjoint.
extern bool FUZZYPARTITION;
///Current name (full path) of partanalyzer
extern char* program;
///Beta parameter for thermodynamic averages done in Ccop
extern double beta;
///mu parameter for thermodynamic averages done in Ccop
extern double mu;
///Original Program basename
#define _programb_ "partanalyzer"

///Type definition of edge
typedef double edge;
///A row is a vector of edges
typedef vector<edge > row;
///A column is the same as a row
typedef row column;
typedef vector< string > svect;
typedef vector< svect > smat;
typedef multiset< string > smset;
typedef set< string > sset;
typedef multimap<string,string > smap;
typedef string node;
typedef map<pair<node,node>,edge > graph ;
typedef pair<string,string> strpair;
typedef pair<string, pair<sset,double> > neighborhood ;
typedef map<string, pair<sset,double> > scover ;
///Type definitions: specific
typedef svect cluster;

///Converts anything to a string
template <typename T>
inline string ToString(const T& x){
        ostringstream os;
        os << x;
        return os.str();
}
///Beta/mu unset/default value allows cheking if user has already set its value or not.
#define BETA_UNSET -666.0
#define MU_UNSET -666.0
#define BETA_DEFAULT 0.01
#define MU_DEFAULT 0.0
///Extensivity degree is precisely that parameter of Tsallis and Renyi entropies. This value means using the default.
#define EXTENSIVITY_DEFAULT -6666.0
#define EXTENSIVITY_DEFAULT_TSALLIS 2.0
#define EXTENSIVITY_DEFAULT_RENYI 2.0
#define EXTENSIVITY_MAX 100.0
//#define RDC_DEFAULT_TOP_BEST_PARTITIONS 20
#define _PGM_P2_GRAYSCALE_ 256
///MSA default format 
#define MSADEFAULTFMT FASTA
///Default Cluster offset
#define CLUSTEROFFSET_DEFAULT 2
///For RDC
#define RDC_DEFAULT_NUMBER_NEIGHBORS 2
///Available types of symmetric and non-symmetric pair-wise measures. 
//enum pmeasure { conditionalEntropy=conditionalPotential , jointEntropy=jointPotential , symmetricPurity , vmeasureArithmetic , vmeasureGeometric , vmeasureHarmonic};
enum pmeasure { conditionalEntropy , jointEntropy , SA, mutualInformation=SA , SSA , conditionalMutualInformation=SSA , SSSA, symmetricPurity , vmeasureArithmetic , vmeasureGeometric , vmeasureHarmonic};
///Available types of pmetric functions (functions inducing a metric on partitions)
enum pmetricv { shannon,entropy=shannon, cardinality , boltzmann, tsallis,renyi, jeffreyQnorm} ;
///Available types of pmetric header comments (#BeginViDistances or #EndEditDistances)
enum flagheader {BEGIN,END};
///Enummerates the possible ways of generating duplicates for a given MSA
enum REDMxVal { useOrgRED, useOwnRED, useZeroRED } ;
///Enummerates possible main program command lines options
enum prganalysis { prgCCOP=1,prgCDIS,
	prgPAXE,prgPASO,prgPASR,prgPACN,
	prgPSPP,prgPSST,prgPSSR,prgPSTG,
	prgVIPP,prgVIST,prgVISR,
	prgEDSC,prgEDST,prgEDSR,
	prgBDST,prgBDSR,
	prgTDST,prgTDSR,
	prgRDST,prgRDSR,
	prgQDST,prgQDSR,
	prgRPAR,prgRPSC,prgRDCB,prgRDCA,prgRPAA,prgRPAB,prgPART,prgPARA,prgPARB,
	prgMEET,prgJOIN,
	prgMEST,prgJOST,
	prgMGMX,prgMMXC, 
	prgCLST,
	prgSNOC,
	prgPMSA,prgPMAS,prgMSPI,prgMAPI,prgMRED,prgMSXP,prgMMAI,prgMMPI,prgMSXS,prgMSDS,prgMSXI,prgMSDI,prgMSXT,prgMSDT,prgMSMP,
	prgHASS,
	prgIPAR,
	prg2MCL,prg2FRE,prgM2PA, prgA2PA,
	prgIPOT,prgIPOR=prgIPOT,prgCPOT,prgJPOT,prgMPOT, prgCKSA=prgMPOT,prgCMPO, prgCSSA=prgCMPO,prgSSSA,prgVMAM,prgVMGM,prgVMHM,prgPSYM,
	prgCPOR,prgJPOR,prgMPOR, prgKSAR=prgMPOR,prgCMPR, prgKSSR=prgCMPR,prgSSSR,prgSEXV,prgVAMR,prgVGMR,prgVHMR,prgPSYR,
	prgIPOG, prgCPOG,prgJPOG,prgMPOG, prgKSAG=prgMPOG,prgCMPG, prgKSSG=prgCMPG,prgSSSG,
	prgADST,prgASFP,
	prgSPST,prgSPSS,prgSPSO,
	prgCEMX,prgEDMX,prgEDMP,prgGPRA,prgGPRB
} ;
///MSA formats: When using --msa-map-partition, both FASTA3 and GDE3 will only output subfamilies with 3 or more sequences.
enum MSAformat { msaFmtNULL=-1, FASTA, GDE , SPEER, FASTA2, GDE2, SPEER2, FASTA3, GDE3, SPEER3, GSIM, GSIM2, GSIM3 } ;
///Enummerates the possible input format for a partition file
enum partFileFormat { partFmtNULL=-1,partFmtPART, partFmtFREE , partFmtMCL, partFmtMSA , partFmtLABEL} ;
///Available types of Splitstat methods
enum splitmethod {split,cosine,overlap,fraction=overlap} ;
///Wrapper of a char* in order to store char* in STL containers. Otherwise, weird memory likages happen, if not simply a segfault.
class Charr{
public: char* car;
};

///Comparison operator for allowing reverse sorting of maps
struct greaterThan {
	template<class T> bool operator() (T& i, T& j){
		return i>j;	
	}
};

template <class C >
struct containerLargerThan {
        bool operator () (C ca, C cb) {
		int offset=1;
		bool condition= ( ca.size()>cb.size() );
		bool equalSize= ( ca.size()==cb.size() );
		if(equalSize)condition=*(ca.begin()+offset) < *(cb.begin()+offset);
		return condition;
        }
};

template <class C>
struct containerLargerThan_Offset {
	int offset;
	//if(offset<0||offset>2)offset=2;
        bool operator () (C ca, C cb) {
		bool condition= ( ca.size()>cb.size() );
		bool equalSize= ( ca.size()==cb.size() );
		if(equalSize)condition=*(ca.begin()+offset) < *(cb.begin()+offset);
		return condition;
        }
};

/**Defines a general operator for beloging (\in) to a set: i.e., does an arbitary element e 
belong to set S? It is assumed that S contains elements of type type(e)=E.
We want to used it as
e < S
Here, we pass the set S by value. Thus we are using a local copy than we can
modify without modifying the original one. We don't need, therefore, to 
remove the element e in case in wasn't there before.
*/
template <class E>
bool operator<(E& e, set<E > S){ 
	long int osize=S.size();
	S.insert(e);
	return (S.size()==osize);
};
	
/**Algebra of string sets
*/
///returns a new set containing the union of sa & sb
template < class E>
set<E > operator+(set<E>& sa, set<E>& sb){
	set<E > suab;
	for(typename set<E>::iterator e=sa.begin();e!=sa.end(); e++)
		suab.insert(*e);
	for(typename set<E>::iterator e=sb.begin();e!=sb.end(); e++)
		suab.insert(*e);
	return suab;
};

///sa+=sb, modifies sa by appending to it sb
template <class E>
set<E>& operator+=(set<E>& sa, set<E>& sb){
	for(typename set<E>::iterator e=sb.begin();e!=sb.end(); e++)
		sa.insert(*e);
	return sa;
	//return (sa=(sa+sb));
};

/* TESTING CUSTOM GRAPH CLASS */
///Definition of sdGraph as a type of graph_base<string, double>
//typedef base_graph<string, double> sd_bGraph;
//
#endif //END OF _PARTANALYZER_DEFINITIONS_HEADER

