#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <exception>
#include <ctime>       /* clock_t, clock, CLOCKS_PER_SEC */
#include <cmath>
#include <algorithm>

#include <boost/graph/topological_sort.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/version.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/concept/assert.hpp>
#include <boost/graph/adjacency_iterator.hpp>
#include <boost/limits.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/undirected_graph.hpp>

using namespace std;
using namespace boost;

struct VertexProperties
{
    int          _index;
    double       _weight;

    VertexProperties() : _index(0), _weight(0)
    { }
    VertexProperties(int index, double weight) :
        _index(index), _weight(weight)
    { }
};

struct EdgeProperties
{
    double       _weight;

    EdgeProperties() : _weight(0)
    { }
    EdgeProperties(double weight) : _weight(weight)
    { }
};

typedef unsigned int uint;
typedef vector<int> Entiers;
typedef vector<Entiers *> EntiersEntiers;
typedef vector<EntiersEntiers> EntiersEntiersEntiers;
typedef list<int> List;
typedef list<EntiersEntiers *> ListEntiersEntiers;
typedef adjacency_list<vecS, vecS, directedS, VertexProperties, EdgeProperties> GraphOriente;
typedef vector < GraphOriente *> Graphs;
typedef adjacency_list<vecS, vecS, undirectedS, VertexProperties, EdgeProperties> GraphNonOriente;
typedef std::pair < int, int > Edge;
typedef vector < Edge > Edges;
typedef Edges OutputEdges;
typedef Edges InputEdges;
typedef vector < OutputEdges > OutputEdgeList;
typedef vector < InputEdges > InputEdgeList;
typedef vector<GraphNonOriente *> Base_Graph;

typedef graph_traits<GraphNonOriente>::vertex_descriptor vertex_t;
typedef graph_traits<GraphNonOriente>::edge_descriptor edge_t;

typedef GraphOriente::vertex_descriptor vertex_to;
typedef GraphOriente::edge_descriptor edge_to;

template <class G> struct VertexAndEdgeListGraphConcept {
	void constraints() {
		function_requires< VertexListGraphConcept<G> >();
		function_requires< EdgeListGraphConcept<G> >();
	}
};


extern GraphNonOriente::vertex_iterator vertexIt, vertexEnd;
extern GraphNonOriente::adjacency_iterator neighbourIt, neighbourEnd;

extern GraphOriente::vertex_iterator vertexIto, vertexEndo;
extern GraphOriente::adjacency_iterator neighbourIto, neighbourEndo;

void Global_Neigh_community(GraphNonOriente *g, const EntiersEntiers &Partition, Entiers *community, int vertex, int comm_in);
Graphs Graph_Partition(const EntiersEntiers &Partition, GraphOriente *go, GraphNonOriente *g, OutputEdgeList &outputedgelist, InputEdgeList &inputedgelist);
void List_edge_partie(Entiers *Partie, GraphOriente *go, Edges &edge_partie, OutputEdges &outputedgespartie);
void construire_graph(GraphNonOriente *g, GraphOriente *graph);
double Modif_Cut_one_cluster(Entiers &cluster, GraphNonOriente &g, double &vol);
vector<double> modif_cut_tmp(GraphNonOriente *g, EntiersEntiers &Partition, vector<vector<double> > tabe_cut, int vertexs, int comm_in, Entiers community, double cut,string name);
double Calcul_poids(Entiers *partie, GraphNonOriente *g);
bool Est_connexe(GraphNonOriente *g, EntiersEntiers Partition, Entiers &part);
void Affinage_equilibrage_charge(GraphNonOriente *g, EntiersEntiers &Partition);
Entiers Neigh_community(GraphNonOriente *g, EntiersEntiers &Partition, int vertex, int comm_in);
void Modif_fonction_Gain_Cut(EntiersEntiers &Partition,GraphNonOriente *g, Entiers &community, int val, double &cut,string name);
void Affinage_recherche_locale(GraphNonOriente *g, EntiersEntiers &Partition, double &cut, string name);
void projection(EntiersEntiers &Partition,ListEntiersEntiers::iterator lit);
void contraction_HEM(GraphNonOriente *g, Base_Graph &baseg, ListEntiersEntiers &liste_corr);
Entiers Liste_adjacence(GraphNonOriente &g, int vertexs, const Entiers &random_vertices);

int rand_fini(int a, int b);
int recherche_val2(const vector<float> &tab,float val);
int recherche_val_double(const vector<double> &tab,double val);
int recherche_val(const vector<int> &tab,int val);
int dichotomie(const Entiers &tab,int i);
void suprim_val(Entiers &tab,int i);
bool In_tab(const Entiers &tab, int val);
bool In_tab_dichotomie(const Entiers &tab, int val);
double Cut_cluster(const EntiersEntiers &tab_cluster,GraphNonOriente &g, string name);
void Modif_fonction_Gain_Cut(EntiersEntiers &Partition,GraphNonOriente *g,int community_out,int community_in,int val,double &cut);
void Liste_Voisin(Entiers &P,Entiers &tab,const GraphNonOriente &g);
int Cout_coupe(Entiers P,int val, GraphNonOriente &g);
double Cout_coupe_pond(Entiers P,int val, GraphNonOriente &g);
int In_community_dichotomie(const EntiersEntiers &part, int val);
double Degree(GraphNonOriente &g , int node);
double Cluster_Degree(GraphNonOriente &g , const Entiers &cluster);

#endif
