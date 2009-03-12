/** Graph 
Copyright (C) Miguel A. Santos, HSC, Toronto, 2008-2009.
Licensed under GPL version 3 a later. (see http://www.gnu.org/copyleft/gpl.html )

Compile Options:


*/

#include<map>
#include<iostream>
/*
#include<vector>
#include<string>
#include<sstream>
#include<set>
#include<list>
#include<math.h>
*/

using namespace std;

/**Abstract base graph
Nodes are of generic type S and edges of type V.
For N nodes, size grows ~ n_pairs * ( 2*sizeOf(int) + sizeOf(V) ) + 2*N*( sizeOf(int) + sizeOf(S) )
thus the dependency on sizeOf(S) is only linearly on N, not with the square of N as in map< pair<S,S> , V>.

*/
template<class S = string, class V = double > class base_graph
{
	map<int,S> _nodes_i;
	map<S,int> _nodes_S;
	map<pair<int, int>, V> _g;
	int _put(S s);
	int _get(S s);
	V& _v(pair<int, int> pi){ return _g[pi]; }
	bool _existEdge(int i, int j);
	pair<S, S>& _cullEdge(int i, int j);
public:
	base_graph();
	base_graph operator= (base_graph&);
	V& operator[] (pair<S, S> s);
	int insert(pair<S,S> s, V v);
	int n_nodes(){return _nodes_S.size();}
	int n_pairs(){return _g.size();}
	bool existEdge(S si, S sj){ return _existEdge(_get(si),_get(sj));}
	bool existNode(S s){ return _get(s)==-1?false:true; }
	pair<S, S> cullEdge(S si, S sj){ return _cullEdge(_get(si),_get(sj));}
	S get(int i);
	V& v(S si, S sj);
	base_graph operator+ (base_graph& ga, base_graph& gb);
	base_graph& operator+= (base_graph& ga, base_graph& gb){ return ga=ga+gb;}
	base_graph operator* (V& z, base_graph& g);
	base_graph operator* (base_graph& g, V& z){ return (z*g) ;}
	ostream& operator<<(ostream& os, base_graph& g);
	class iterator;
	class iterator_index;
	class iterator_node;
};

template<class S, class V> base_graph<S,V>::base_graph(){
}

template<class S, class V> int base_graph<S,V>::_put(S s){
	int n=_nodes_S.size();
	_nodes_S[s]=n;
	_nodes_i[n]=s;
	return n;
}

template<class S, class V> int base_graph<S,V>::_get(S s){
	if(_nodes_S.find(s)!=_nodes_S.end()) return _nodes_S[s];
	return -1;
}

template<class S, class V> base_graph& base_graph<S,V>::operator=(base_graph g){
	_nodes_i=g._nodes_i;
	_nodes
}

template<class S, class V> V& base_graph<S,V>::operator[](pair<S, S> s){
	if(!existEdge(s.first,s.second)) insert(s, (V)NULL);
	return _v( pair<int, int> (_get(s.first),_get(s.second)) ) ;
}

template<class S, class V> V& base_graph<S,V>::v(S si, S sj){
	if(existEdge(si,sj)) return _v(pair<int, int> (get(si),get(sj))) ;
}

template<class S, class V> S base_graph<S,V>::get(int i){
	if(_nodes_i.find(i)!=_nodes_i.end()) return _nodes_i[i];
	return (S) NULL;
}

template<class S, class V> bool base_graph<S,V>::_existEdge(int i, int j){
	if(i<0||j<0) return false;
	if(_g.find(pair<int, int> (i,j))!=_g.end()) return true;
	return false;
}

template<class S, class V> pair<S, S>& base_graph<S,V>::_cullEdge(int i, int j){
	if(_existEdge(i,j)) return _g[pair<int, int> (i,j)];
	return  &pair<S, S> (NULL,NULL);
}

template<class S, class V> int base_graph<S,V>::insert(pair<S, S> s, V v){
	if(_nodes_S.find(s.first)==_nodes_S.end()) _put(s.first);
	if(_nodes_S.find(s.second)==_nodes_S.end()) _put(s.second);
	pair<int, int> pi (_get(s.first),_get(s.second));
	_g[pi]=v;
	return n_pairs();
}

template<class S, class V> base_graph base_graph<S,V>::operator+(base_graph ga, base_graph gb){
	base_graph sg;
	base_graph lg=ga;
	if(ga.n_nodes()<gb.n_nodes()) lg=gb;
	for(int i=0;i<lg;i++)
		
	return sg;
}

class base_graph<S,V>::iterator{

};

class base_graph<S,V>::iterator_index{

};

class base_graph<S,V>::iterator_node{

};

