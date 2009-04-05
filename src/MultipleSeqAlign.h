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

#ifndef _CLASS_MULTIPLESEQALIGN_H
#define _CLASS_MULTIPLESEQALIGN_H 1

#include "Sequence.h"
#include "Roulette.h"

///Type definition of a MSA as a vector of Sequences
typedef vector< Sequence > MSA;
///**Allows to perform operations on multiple sequence aligments. Its main member being of course the MSA object _Seqlist which
//contains the list of all sequences and their names.
//*/
class MultipleSeqAlign
{
        char* _msaf;
	string _name; ///General purpose label
        MSA _Seqlist;
        int _len;
        int _nseq;
	roulette ruleta;
	//double _id_thr; ///For calculating Id statistics: id threshold.

public:
	///Default constructor
        MultipleSeqAlign();
	///Expects a MSA file in fasta format
        MultipleSeqAlign(char* msaf); 
	///Expects a MSA object
	MultipleSeqAlign(MSA& msa);
	///Subsample MSA at given positions
	MultipleSeqAlign xtractPositions(vector<int>* positions=NULL);
	///Extract sequences specified by name. IF bool equal=false, drop them instead.
	MultipleSeqAlign xtractSequences(svect* seqnames, bool equal=true);
	///Get begin iterator to _Seqlist
	MSA::iterator beginSeq(){ return _Seqlist.begin();}
	///Get end iterator to _Seqlist
	MSA::iterator endSeq(){ return _Seqlist.end();}
	///Allows adding a single Sequence object to the multiple sequence alignment 
	void addSeq(Sequence*);
	void addSeq(Sequence);
	///Is MSA empty?
	bool empty(){ return _Seqlist.empty(); }
	///Sets its name
	void setName(string name){ _name=name;}
	///Gets its name
	string getName(){ return _name;}
	/**Sets its file name. This, being just a char* can be used as an additional lable. That can be helpful with multiple sequence alignments
	generated on the fly, during the calculations.*/
	void setFileName(char* fname){ _msaf=fname;}
	///Gets its file name
	char* getFileName(){ return _msaf;}
	///Obtains the sequence identity between the two provided Sequences using the specified positions.
        double SeqId(Sequence Seqa, Sequence Seqb, vector<int>* positions=NULL){return Seqa.id(Seqb,positions);}
        //double SeqId(string SeqaName, string SeqbName);
	///Obtains the sequence identity between the two provided sequences refered by their indexes within the MSA.
        double SeqId(int Seqn, int Seqm,vector<int>* positions=NULL);
	///Calculates the overall average sequence identity among all pair of sequences of the MSA.
        double averageId(vector<int>* positions=NULL);
	///Prints all pair-wise identities 
        void printPairwiseIds(vector<int>* positions=NULL);
	///Prints the average pair-wise identity and the fraction of all pairs with identity above thr. Default thr=50
        void printAveragePairwiseIds(double thr=50.0,vector<int>* positions=NULL);
	///Given alternative MSA, prints all pair-wise identities 
        void printPairwiseIds(MultipleSeqAlign& msa, vector<int>* positions=NULL);
	///Given alternative MSA, prints the average pair-wise identity and the fraction of all pairs with identity above thr. Default thr=50. 
        void printAveragePairwiseIds(MultipleSeqAlign& msa, double thr,vector<int>* positions=NULL);
	///Print the whole multiple sequence alignment
        void print();
	///Generates randomly n additional sequences each an exact copy of one of the original sequences
	MultipleSeqAlign genRedundantMSA(int n);
	/**Generaes nsamples of MSAs each with redundancy of size n, i.e., each having n redundant sequences. The system clock is used for 
	generating the seed of the random number generator*/
	void NsamplesRedMSA(int nsamples, int n);
	///The same as the previous one, but allows specifying a particular seed
	void NsamplesRedMSA(unsigned int seed,int nsamples, int n);
	//Streams out a multiple sequence alignments (in fasta format) wihtout adding a new line at the end
	friend ostream& operator<<(ostream& os, MultipleSeqAlign& msa);
};

#endif //END _CLASS_MULTIPLESEQALIGN_H

