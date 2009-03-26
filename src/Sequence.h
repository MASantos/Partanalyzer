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

/**Allows dealing with sequences (protein sequences, DNA sequences,etc.).
 * */
class Sequence
{
        string _seq;
        string _name;
        int _length;
public:
	///Default constructor
        Sequence();
	///Instatiate a sequence given its name label and the sequence proper seq
        Sequence(string name, string seq);
	///Get the character at position pos
        string chrAt(int pos);
	///Get its lenght
        int length(){ return _length;}
	///Calculate its sequence identity with respect to sequence Seq. The latter can be passed by pointer or content
        ///A non null second parameter specifies which positions will be considered to calculate the sequence similarity, all by default.
        double id(Sequence* Seq, vector<int>* positions=NULL);
        double id(Sequence Seq, vector<int>* positions=NULL);
	///Get its name
        string name(){return _name;}
	///Get its sequence proper
        string sequence(){return _seq;}
	///Set its sequence 
        void setSeq(string seq){ _seq=seq; _length=_seq.length();}
	///Set its name
        void setName(string name){ _name=name;}
	///Print a sequence with its name adding a newline at the end
        void print(){cout<<_name<<"\n"<<_seq<<endl;}
	///Streams out a sequence and its name. Does not add a newline at the end.
	friend ostream& operator<<(ostream& os , Sequence& s);
};
#endif //END _CLASS_SEQUENCE_H

