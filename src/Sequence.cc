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

#ifndef _CLASS_SEQUENCE
#define _CLASS_SEQUENCE 1

#include "Sequence.h"

pid_normalization PIDNORMALIZATION=shorterSequence;

void printPIDNormalization(){
	cout<<"#PIDNORMALIZATION= ";
	switch(PIDNORMALIZATION){
		case shorterSequence:
			cout<<"shorterSequence";
			break;
		case arithmeticMeanSequenceLenth:
			cout<<"arithmeticMeanSequenceLenth";
			break;
		case numberOfAlignedResiduePairs:
			cout<<"numberOfAlignedResiduePairs:";
			break;
		case numberOfAlignedPositions://alignmentLenth:
			cout<<"numberOfAlignedPositions";
			break;
	}
	cout<<endl;
}

ostream& operator<<(ostream& os , Sequence& s){ os<<s._name<<"\n"<<s._alignedsequence; return os;}

Sequence::Sequence()
{
        _name="";
        _alignedsequence=_sequence="";
        _alignment_length=_number_of_residues=0;
	_first_residue_number_gauge=0;
}

Sequence::Sequence(string name, string seq){
        _alignedsequence=seq;
        _name=name;
        _alignment_length=_alignedsequence.length();
	_setActualSeq();
}

void Sequence::_setActualSeq(){
	removeGaps();	
	_number_of_residues=_sequence.length();
}

void Sequence::setAlignedSeq(string seq){
	 _alignedsequence=seq; 
	_alignment_length=_alignedsequence.length();
	_setActualSeq();
}

string Sequence::removeGaps(){
	string seqwogaps;
	for(int i=0;i<_alignedsequence.length();i++){
		if(_alignedsequence.substr(i,1).compare("-")==0)continue;
		seqwogaps+=_alignedsequence.at(i);
	}
	_sequence=seqwogaps;
	return seqwogaps;
}

int Sequence::firstResidueNumber(int nb){
	if(nb=0){
		return _first_residue_number_gauge;
	}
	else{
		_first_residue_number_gauge=nb-1;
		return _first_residue_number_gauge+1;
	}
}

string Sequence::chrAt(int pos){
        if(pos<0||pos>_alignment_length-1){
                cout<<"Seeking character beyond sequence boundaries "<<_alignment_length<<" : "<<pos<<endl;
                exit(1);
        }
        return _alignedsequence.substr(pos,1);
}

string Sequence::xtractPositions(vector<int>* positions){
	string nseq="NAS";
	if(positions->empty())return nseq;
	stringstream ss;
	for(vector<int>::iterator it=positions->begin();it!=positions->end();it++){
		for(int i=0;i<_alignment_length;i++){
			if(i==*it) ss<<chrAt(i);
		}
	}
	ss>>nseq;
	return nseq;
}

double Sequence::_getPIDNormalization(Sequence* Seq, int& starts, int& ends){
	switch(PIDNORMALIZATION){
		case shorterSequence:{
			int minlen=Seq->numberOfResidues();
			if(numberOfResidues()<minlen)minlen=numberOfResidues();
			return minlen*1.0;
			break;
			}
		case arithmeticMeanSequenceLenth:
			return 0.5*(numberOfResidues()+Seq->numberOfResidues());
			break;
		case numberOfAlignedResiduePairs:
		case numberOfAlignedPositions:{//alignmentLenth:
			int nb_alignedpos=0;
			int nb_alignedres=0;
			starts=0; //Points to one position BEFORE the left most residue of either sequence (paired or unpaired)
			ends=_alignment_length;//Points to one position AFTER the right most residue of either sequence (paired or unpaired)
        		bool founds=false;
        		bool founde=false;
        		for(int i=0;i<_alignment_length;i++){
        		        if(!founds&& (chrAt(i).compare("-")==0&&Seq->chrAt(i).compare("-")==0) )starts++;
        		        else founds=true;
        		        if(!founde&& (chrAt(_alignment_length-1-i).compare("-")==0&&Seq->chrAt(_alignment_length-1-i).compare("-")==0) )ends--;
        		        else founde=true;
        		        if(!founds)continue; //Calculate nb_alignedpos and nb_alignedres once found the left most residue
                		if(chrAt(i).compare("-")==0&&Seq->chrAt(i).compare("-")==0)continue;
				nb_alignedpos++;
				//If the position does not contain a gap in either sequence, it is an aligned residue pair
				if(!(chrAt(i).compare("-")==0||Seq->chrAt(i).compare("-")==0))
					nb_alignedres++;
        		}
			if(VERBOSE)cout<<"#numberOfAlignedPositions= "<<nb_alignedpos<<" numberOfAlignedResiduePairs= "<<nb_alignedres<<endl;
			switch(PIDNORMALIZATION){
				case numberOfAlignedResiduePairs:
					return nb_alignedres*1.0;
					break;
				case numberOfAlignedPositions:
					return nb_alignedpos*1.0;
					break;
			}
			break;
			}
	}
}
double Sequence::id(Sequence Seq, vector<int >* positions){
        id(&Seq, positions);
}
double Sequence::id(Sequence* Seq,vector<int>* positions)
{
        int other_len=Seq->alignmentLength();
        int minlen=_alignment_length<other_len?_alignment_length:other_len;
        if(minlen!=_alignment_length){
                cout<<"Comparing sequence alignments of different lengths: "<<_name<<" ("<<_alignment_length<<") <> "<<Seq->name()<<" ("<<other_len<<")"<<endl;
                exit(1);
        }
        int starts=0;//Points to one position BEFORE the left most residue of either sequence (paired or unpaired)
        int ends=_alignment_length;//Points to one position AFTER the right most residue of either sequence (paired or unpaired)
        bool all_positions=true;
	//Sequence ID normalization
        double norm=_getPIDNormalization(Seq,starts,ends);
        if(DEBUG)cout<<"#Sequence comparison starting at position "<<starts+1<<" up to position "<<ends<<endl;
        int selected_pos=-9999;
        double id=0.0;
        for(int i=starts;i<ends;i++){
#ifdef DEBUG
		if(VERBOSE)cout<<" starts="<<starts<<" ends="<<ends<<" : i="<<i<<endl;
#endif
                if(!(positions==NULL || positions->empty() ) ){
                        all_positions=false;
                        for(vector<int>::iterator c=positions->begin();c!=positions->end();c++){
                                if(i==*c){
                                        selected_pos=i;
                                        break;
                                }
                        }
                        if(selected_pos!=i)continue;
                        if(DEBUG)cout<<"#Position: "<<i<<endl;
                }
                if(chrAt(i).compare("-")==0&&Seq->chrAt(i).compare("-")==0)continue;
                if(chrAt(i).compare(Seq->chrAt(i))==0)id++;
        }
        //return id/(ends-starts)*100.0;
        return id/norm*100.0;
}
/*
double Sequence::id(Sequence Seq){
        id(&Seq);
}
double Sequence::id(Sequence* Seq)
{
        int other_len=Seq->length();
        int minlen=_length<other_len?_length:other_len;
        if(minlen!=_length){
                cout<<"Comparing sequence of different lengths: "<<_name<<" ("<<_length<<") <> "<<Seq->name()<<" ("<<other_len<<")"<<endl;
                exit(1);
        }
        double id=0.0;
        int starts=0;
        int ends=_length;
        bool founds=false;
        bool founde=false;
        for(int i=0;i<_length;i++){
                /// Normalize by the longest sequence
                /<
                if(!founds&&chrAt(i).compare("-")==0&&Seq->chrAt(i).compare("-")==0)starts++;
                else founds=true;
                if(!founde&&chrAt(_length-1-i).compare("-")==0&&Seq->chrAt(_length-1-i).compare("-")==0)ends--;
                else founde=true;
                >/
                /// Normalize by the shortest sequence
                if(!founds&& (chrAt(i).compare("-")==0||Seq->chrAt(i).compare("-")==0) )starts++;
                else founds=true;
                if(!founde&& (chrAt(_length-1-i).compare("-")==0||Seq->chrAt(_length-1-i).compare("-")==0) )ends--;
                else founde=true;
                //
                if(founds&&founde)break;

        }
        double norm=0;
        if(VERBOSE)cout<<"#Sequence comparison starting at position "<<starts+1<<" up to position "<<ends<<endl;
        for(int i=starts;i<ends;i++){
                if(chrAt(i).compare("-")==0&&Seq->chrAt(i).compare("-")==0)continue;
                norm++;
                if(chrAt(i).compare(Seq->chrAt(i))==0)id++;
                //if(strcmp(chrAt(i).c_str(),Seq->chrAt(i).c_str())==0)id+=1.0;
        }
        //return id/_length*100.0;
        //return id/(ends-starts)*100.0;
        return id/norm*100.0;
}
*/

void Sequence::printAlignment(bool withoutgaps){
	if(_name.compare("")==0||_alignedsequence.compare("")==0)return;
	if(withoutgaps){
		cout<<_name<<"\n"<<_sequence<<endl;
	}else{
		cout<<_name<<"\n"<<_alignedsequence<<endl;
	}
}

#endif //END _CLASS_SEQUENCE

