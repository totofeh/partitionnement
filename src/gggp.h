#ifndef GGGP_H
#define GGGP_H

#include "util.h"

void gggp(GraphNonOriente *g,Entiers *sommetsSource, Entiers *sommetsDestination,EntiersEntiers &Partition);
void gggp_pond(GraphNonOriente *g,Entiers *sommetsSource, Entiers *sommetsDestination,EntiersEntiers &Partition);
void Iter_2l(EntiersEntiers &part, int nbr_parties, GraphNonOriente *g,const string &nom);
void bissectionRec(GraphNonOriente *g,EntiersEntiers &Partition,int nbr_parties,const string &nom);
void Pseudo_random_partitioning(GraphNonOriente *g, EntiersEntiers &Partition, int nbr_parties);
Graphs Multiniveau(int niveau_contraction, GraphNonOriente *g, GraphOriente *go, int nbr_parties, string type_methode, string choix_affinage, string type_cut,Edges &edge_partie , OutputEdgeList &outputedgeslist, InputEdgeList &inputedgelist, Connections &connections);

#endif

