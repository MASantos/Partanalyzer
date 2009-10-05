#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
#include <string>
#include <dirent.h>
#include <sys/stat.h>
#include <vector>
using namespace std;

string getSubstr(string a){
	a=a.substr(0,a.find("-RED"));
	return a;
}

/*
template<class T> class kk {
public:
	map<int, T> k;
	T& operator[] (int i){ return k[i];}
	class iterator;
	iterator& begin(){ return k.begin() ;}
	iterator& end(){ return k.end() ;}
};

template<class T> class kk::iterator : public std::iterator<std::random_access_iterator_tag,kk>{
	map<int, T>& mymap;
public:
	iterator():mymap(){}
	iterator(const map<int, T>& mup):mymap(mup){}
	
};
*/

//bool belongsto(E& e, set<E > S){
template <class E>
bool operator<(E& e, set<E > S){
        long int osize=S.size();
	S.insert(e);
        return (S.size()==osize);
};

///returns a new set containing the union of sa & sb
template < class E>
set<E> operator+(set<E>& sa, set<E>& sb){
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
        return (sa=(sa+sb));
};

int main(){
	int intarray[] = {1,2,3};
	set<int > s(intarray,intarray+3);
	int ia=4;
	cout<<"set s has "<<s.size()<<" members"<<endl;
	//if(!belongsto(ia,s))cout<<ia<<" does not belong to s"<<endl;
	if(!(ia<s))cout<<ia<<" does not belong to s"<<endl;
	else cout<<ia<<" belongs to s"<<endl;
	cout<<"set s has "<<s.size()<<" members"<<endl;
	int intarray2[] = {4,5};
	set<int > s2(intarray2,intarray2+2);
	cout<<"set s2 has "<<s2.size()<<" members"<<endl;
	set<int > ss=s+s2;
	cout<<"set s+s2 has "<<ss.size()<<" members"<<endl;
	/*
	map<string , double > mymap;
	mymap.insert(pair<string,double> ("a",0.988));
	mymap.insert(pair<string,double> ("b",0.0));
	string sa,sb;
	sa="a";
	sb="b";
	cout<<"printing:"<<endl; 
	cout<<sa<<" "<<mymap[sa]<<endl;
	cout<<sb<<" "<<mymap[sb]<<endl;
	if(mymap[sb])cout<<sb<<" "<<"Eureka"<<endl;
	else cout<<"Oops..."<<endl;
	sb="c";
	cout<<sb<<" "<<mymap[sb]<<endl;
	if(mymap[sb])cout<<sb<<" "<<"Eureka"<<endl;
	else cout<<"Oops..."<<endl;
	double dd=1.432;
	double de=1.532;
	cout<<"dd="<<dd<<" de="<<de<<endl;
	int ide=(int)de;
	int idd=(int)dd;
	cout<<"idd="<<idd<<" ide="<<ide<<endl;
	//mymap.insert(pair<string,double> (sb,de));
	mymap[sb]=de;
	cout<<"#map after adding pair ("<<sb<<","<<de<<")"<<endl;
	for(map<string , double >::iterator it=mymap.begin();it!=mymap.end();it++)
		cout<<"("<<it->first<<","<<it->second<<")"<<endl;

	cout<<endl;
	int er=system(NULL);
	cout<<"SHELL exists: shell returned with value="<<er<<endl;
	int odir=system("mkdir cagada");
	if(odir!=0){
		cout<<"shell returned with value="<<er<<endl;
		cout<<"Failed to open dir"<<endl;
	}
	else {
		cout<<"shell returned with value="<<er<<endl;
		cout<<"directory opened ok"<<endl;
	}
	ofstream of("cagada/two");
	if(!of)cout<<"Failed to open file cagada/two"<<endl;
	else of<<"cagada/two"<<endl;
	of.close();
	system("ls cagada");
	*/
	//char* pdir="/home/shoshana/msantos/DNA-PROTEIN/PRG/wodaks/partanalyze/casta";
	


	/*
	char* pdir="casta";
	//unsigned int mode=0766;
	mode_t mode=0744;
	int od=mkdir(pdir,mode);
	if(od!=0)cout<<"Failed to create directory "<<pdir<<endl;
	string fn="one";
	stringstream ss;
	ss<<"./"<<pdir<<"/"<<fn;
	ofstream oft(ss.str().c_str());
	//of(ss.str().c_str());
	if(!oft)cout<<"Failed to open file "<<ss.str()<<endl;
	else oft<<"ok"<<endl;
	*/

	/*string a("Mi mama-RED");
	cout<<"Initially a="<<a<<endl;
	cout<<"function returns="<<getSubstr(a)<<endl;
	cout<<"We still see a="<<a<<endl;
	ifstream is("alpha.txt");
	string st,fc;
	int isint;
	getline(is,st,'$');
	cout<<"The first line seen is: ->"<<st<<"<-"<<endl;
	stringstream ss;
	ss<<st;
	cout<<"printing its words..."<<endl;
	while(ss>>st )if(st.compare("\n")==0){cout<<"EOL"<<endl;}else{cout<<st<<endl;}
	cout<<"Done with that line."<<endl;
	//while(is>>st && st.substr(0,1).compare("#")!=0 && isint=atoi( &st.at(0) ) && 0<=isint && isint<=9 )
	while(is>>st )
		if(st.substr(0,1).compare("#")!=0 && isdigit(st.at(0)) ) 
			cout<<st<<" starts with a digit"<<endl;
		//else
		//	cout<<st<<" doesn't start with a digit"<<endl;
	*/
	/*
	*/
	vector< int > v;
	cout<<"Declaring vector: v={";
	for(int i=0;i<10;i++){
		v.push_back(i+1);
		cout<<v[i]<<",";
	}
	cout<<"}"<<endl;
	vector< vector< int > > m;
	m.push_back(v);
	cout<<"erasing 6th element"<<endl;
	vector< int>::iterator itv=v.begin();
	v.erase(itv+5);
	//As iterator
	vector< vector< int> >::iterator vit=m.begin();
	vector< int>::iterator it=vit->begin();
	cout<<"erasing element #1: "<<*it<<endl;
	vit->erase(it);
	cout<<"v={";
	for(int i=0;i<v.size();i++)cout<<v[i]<<",";
	cout<<"}"<<endl;
	cout<<"m[0]={";
	for(int i=0;i<m[0].size();i++)cout<<m[0][i]<<",";
	cout<<"}"<<endl;
	string stra="I am 21 characters long.";
	string strb="Four";
	string strc="I am 84 characters long.I am 84 characters long.I am 84 characters long.I am 84 characters long.";
	map<pair<string,string>,double> mymap;
	mymap[pair<string,string> (stra,stra)]=1.0;
	mymap[pair<string,string> (stra,strb)]=2.0;
	mymap[pair<string,string> (strc,strb)]=3.0;
	cout<<"stra="<<stra<<endl;
	cout<<"strb="<<strb<<endl;
	cout<<"strc="<<strc<<endl;
	cout<<"sizeof(stra)="<<sizeof(stra)<<endl;
	cout<<"sizeof(strb)="<<sizeof(stra)<<endl;
	cout<<"sizeof(strc)="<<sizeof(stra)<<endl;
	cout<<"sizeof(int)="<<sizeof(int)<<endl;
	cout<<"sizeof(long int)="<<sizeof(long int)<<endl;
	cout<<"sizeof(double)="<<sizeof(double)<<endl;
	cout<<"sizeof(*str)="<<sizeof(&stra)<<endl;
	cout<<"sizeof(mymap)="<<sizeof(mymap)<<endl;
	/*
	cout<<"Instantiating template class kk"<<endl;
	kk<int> mmap;
	mmap[0]=1;
	cout<<"First element : mmap[0]= "<<mmap[0]<<endl;
	*/
	return 0;	
}
