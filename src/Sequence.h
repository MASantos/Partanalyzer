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

#ifndef _CLASS_SEQUENCE_H
#define _CLASS_SEQUENCE_H 1

#include "partanalyzer_includes.h"
#include "partanalyzer_definitions.h"

/**See Alex C.W. May, "Percent Sequenc Identity: The need to be explicit", Structure, vol. 12, 737-738, May, 2004
numberOfalignedPositions: includes gaps
numberOfAlignedResiduePairs: does not include gaps
*/
enum pid_normalization { 
	shorterSequence,
	alignmentLenth,
	numberOfAlignedPositions=alignmentLenth,
	numberOfAlignedResiduePairs,
	arithmeticMeanSequenceLenth
};

extern pid_normalization PIDNORMALIZATION;
extern void printPIDNormalization();

/**Allows dealing with sequences (protein sequences, DNA sequences,etc.).
 * */
class Sequence
{
	///Name of sequence
        string _name;
	///String of characters, including gaps. Corresponds the string as in a MSA.
        string _alignedsequence;
	///Format of aligned sequence (MSA format)
	MSAformat _msafmt;
	///Number of characters, including gaps.  Corresponds the length of a MSA.
        int _alignment_length;
	///Actual string of characters, without gaps
        string _sequence;
	///Actual Number of characters, without gaps
        int _number_of_residues;
	///Sets actual sequence string and actual length base on _alignedsequence
	void _setActualSeq();
	///Residue number offset: gauges the actual residue numbers
	int _first_residue_number_gauge;
	///Get the PID normalization for the given two sequences
	double _getPIDNormalization(Sequence* Seq, int& starts, int& ends);
	///Set sequence's _msafmt variable depending on the format of its _name
	void _setFormat();
public:
	///Default constructor
        Sequence();
	///Instatiate a sequence given its name label and the sequence proper seq
        Sequence(string name, string seq);
	///Sets actual residue number of first residue and return it. If nb=0, just returns first residue number 
	int firstResidueNumber(int nb=0);
	///Get the character at position pos of the MSA sequence
        string chrAt(int pos);
	///Get the character at position pos of the actual sequence
        string chrAt_bare(int pos);
	///Copy extracting from sequence the given positions
	string xtractPositions(vector<int>* positions=NULL);
	///Get its MSA lenght
        int alignmentLength(){ return _alignment_length;}
	///Get its actual lenght
        int numberOfResidues(){ return _number_of_residues;}
	///Calculate its sequence identity with respect to sequence Seq. The latter can be passed by pointer or content
        ///A non null second parameter specifies which positions will be considered to calculate the sequence similarity, all by default.
        double id(Sequence* Seq, vector<int>* positions=NULL);
        double id(Sequence Seq, vector<int>* positions=NULL);
	///Get its name
        string name(){return _name;}
	///Get its sequence as in MSA, including possible gaps
        string alignedSequence(){return _alignedsequence;}
	///Get its actual sequence, without gaps
        string sequence(){return _sequence;}
	///Set its sequence 
        void setAlignedSeq(string seq);
	///Set its name
        void setName(string name){ _name=name;}
	///Stream out its name in the given MSA format. Does not actually change the sequence's format
	string formatName(MSAformat msafmt=MSADEFAULTFMT);
	///Set MSA format. This does change the sequence's format
	void setFormat(MSAformat msafmt);
	///Print a sequence with its name adding a newline at the end. If bare=true, print its actual sequence, i.e., without gaps
        void printAlignment(bool withoutgaps=false);
        void printAlignment(MSAformat msafmt, bool withoutgaps=false);
	///Print the actual sequence with its name adding a newline at the end
        void print(){ return printAlignment(true);}
	///Remove gaps 
	string removeGaps();
	///Maps actual residue numbers to aligned position numbers. If NULL, returns MSA positions of all residues (gaps ignored)
	vector<int > getMSAPositionsFromActualResidueNumbers(vector<int>* residueNumbers=NULL);
	///Maps aligned position numbers to actual residue numbers. If NULL, returns Actual residue numbers of all residues.
	vector<int > getActualResidueNumbersFromMSAPositions(vector<int>* residueNumbers=NULL);
	///Streams out a sequence and its name. Does not add a newline at the end.
	friend ostream& operator<<(ostream& os , Sequence& s);
};
#endif //END _CLASS_SEQUENCE_H

