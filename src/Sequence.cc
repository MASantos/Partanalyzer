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

ostream& operator<<(ostream& os , Sequence& s){ os<<s._name<<"\n"<<s._seq; return os;}

Sequence::Sequence()
{
        _seq="";
        _name="";
        _length=0;
}
Sequence::Sequence(string name, string seq){
        _seq=seq;
        _name=name;
        _length=_seq.length();
}
string Sequence::chrAt(int pos){
        if(pos<0||pos>_length-1){
                cout<<"Seeking character beyond sequence boundaries "<<_length<<" : "<<pos<<endl;
                exit(1);
        }
        return _seq.substr(pos,1);
}

double Sequence::id(Sequence Seq, vector<int >* positions){
        id(&Seq, positions);
}
double Sequence::id(Sequence* Seq,vector<int>* positions)
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
        bool all_positions=true;
        for(int i=0;i<_length;i++){
                /// Normalize by the longest sequence
                /*
                if(!founds&&chrAt(i).compare("-")==0&&Seq->chrAt(i).compare("-")==0)starts++;
                else founds=true;
                if(!founde&&chrAt(_length-1-i).compare("-")==0&&Seq->chrAt(_length-1-i).compare("-")==0)ends--;
                else founde=true;
                */
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
        int selected_pos=-9999;
        for(int i=starts;i<ends;i++){
                //if(!(positions==NULL)&&positions->empty()>0){
                if(!positions->empty()){
                        all_positions=false;
                        for(vector<int>::iterator c=positions->begin();c!=positions->end();c++){
                                if(i==*c){
                                        selected_pos=i;
                                        break;
                                }
                        }
                        if(selected_pos!=i)continue;
                        if(VERBOSE)cout<<"#Position: "<<i<<endl;
                }else exit(1);
                if(chrAt(i).compare("-")==0&&Seq->chrAt(i).compare("-")==0)continue;
                norm++;
                if(chrAt(i).compare(Seq->chrAt(i))==0)id++;
                //if(strcmp(chrAt(i).c_str(),Seq->chrAt(i).c_str())==0)id+=1.0;
        }
        //return id/_length*100.0;
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
#endif //END _CLASS_SEQUENCE

