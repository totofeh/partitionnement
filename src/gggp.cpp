#include "gggp.h"

using namespace std;
using namespace boost;

void gggp(GraphNonOriente *g,Entiers *sommetsSource, Entiers *sommetsDestination,EntiersEntiers &Partition)
{
	int val;
	Entiers sommets_adj;
	if((sommetsSource->size()-1)==0)
	{
		val=0;
		//cout<<"Entré dans le debug ! "<<endl;
		Entiers tailles;

		for(int i=0;i<Partition.size();i++)
		{
			tailles.push_back(Partition.at(i)->size());
		}
		int tmp=*max_element(tailles.begin(),tailles.end());
		for(int i=0; i<Partition.size();i++)
		{
			if(Partition.at(i)->size()==tmp)
				gggp(g,Partition.at(i),sommetsDestination,Partition);
			break;
		}
	}
	else
		val=rand_fini(0,sommetsSource->size()-1);//Tirage aléatoire de l'indice du premier sommet entre 0 et taille du tableau -1

	float poids_max=sommetsSource->size()/2.;
	float poids=1;
	Entiers sommets_cut;

	//clog<<"Etape 1 : "<<endl;
	sommetsDestination->push_back(sommetsSource->at(val));
	sommetsSource->erase(sommetsSource->begin() + val);

	if(sommetsSource->size()<2)
		return;

	while(poids<poids_max)
	{
//		for(uint i =0; i< sommetsDestination.size();i++){
//			cout<<sommetsDestination.at(i)<<endl;
//		}
		Liste_Voisin(*sommetsDestination,sommets_adj,*g);
		if((sommets_adj.size()==0))
		{
			cout<<"Je suis sorti !!!! "<<endl;
			break;
		}
		else{
			/*clog<<"Liste voisin est : "<<endl;
			for(int i=0;i<sommets_adj.size();i++)
			{
				cout<<sommets_adj[i]<<endl;
			}*/
			sort(sommets_adj);
			for(int i=0;i<sommets_adj.size();i++)
			{
				sommets_cut.push_back(Cout_coupe(*sommetsDestination,sommets_adj[i],*g));
			}
			int tmp = recherche_val(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()));
			sommetsDestination->push_back(sommets_adj[tmp]);
			suprim_val(*sommetsSource, sommets_adj[tmp]);
			suprim_val(sommets_adj, sommets_adj[tmp]);

			sommets_cut.clear();
			poids++;

		}
	}

	for (int i=0; i<sommetsSource->size();i++)
	{
		for (int j=0; j<sommetsDestination->size();j++)
		{
			remove_edge(sommetsSource->at(i),sommetsDestination->at(j),*g);
		}
	}

	sort(*sommetsDestination);
}

void gggp_pond(GraphNonOriente *g,Entiers *sommetsSource, Entiers *sommetsDestination,EntiersEntiers &Partition){
	int val;
	Entiers sommets_adj;
	if((sommetsSource->size()-1)==0){
		val=0;
		cout<<"Entré dans le debug ! "<<endl;
		Entiers tailles;
		for(int i=0;i<Partition.size();i++){
			tailles.push_back(Partition.at(i)->size());
		}
		int tmp=*max_element(tailles.begin(),tailles.end());
		for(int i=0; i<Partition.size();i++){
			if(Partition.at(i)->size()==tmp)
				gggp_pond(g,Partition[i],sommetsDestination,Partition);
			break;
		}
	}
	else
		val=rand_fini(0,sommetsSource->size()-1);//Tirage aléatoire de l'indice du premier sommet entre 0 et taille du tableau -1
	double poids_max=0;
	for(int i=0;i<sommetsSource->size();i++){
		poids_max+=(*g)[sommetsSource->at(i)]._weight;
	}
	poids_max/=2.;
	double poids=(*g)[sommetsSource->at(val)]._weight;
	vector<float> sommets_cut;

	sommetsDestination->push_back(sommetsSource->at(val));
	sommetsSource->erase(sommetsSource->begin() + val);

	if(sommetsSource->size()<2)
			return;

	while(poids<poids_max)
	{
		Liste_Voisin(*sommetsDestination,sommets_adj,*g);
		if((sommets_adj.size()==0))
		{
			cout<<"Je suis sorti !!!! "<<endl;
			return;
		}
		else{
			sort(sommets_adj);

			for(int i=0;i<sommets_adj.size();i++)
			{
				sommets_cut.push_back(Cout_coupe_pond(*sommetsDestination,sommets_adj[i],*g));
			}
			sommetsDestination->push_back(sommets_adj[recherche_val2(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]);
			poids+=(*g)[sommets_adj[recherche_val2(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]]._weight;
			suprim_val(*sommetsSource, sommets_adj[recherche_val2(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]);
			suprim_val(sommets_adj, sommets_adj[recherche_val2(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]);

			sommets_cut.clear();

		}
	}

	edge_t e1;
	bool found;

	for (int i=0; i<sommetsSource->size();i++)
	{
		for (int j=0; j<sommetsDestination->size();j++)
		{
			//tie(e1,found)=edge(vertex(sommetsSource->at(i),*g),vertex(sommetsDestination->at(j),*g),*g);
			//if((*g)[e1]._weight!=0){
				remove_edge(sommetsSource->at(i),sommetsDestination->at(j),*g);
			//}
		}
	}

	sort(*sommetsDestination);
}

void Iter_2l(EntiersEntiers &part, int nbr_parties, GraphNonOriente *g,const string &nom)
{
	if (nom!="gggp_pond"){
		//cout<<"je jsuis dans gggp"<<endl;

		for(int i = 0; i<floor(log(nbr_parties)/log(2)); i++)
		{
			//cout<<"Et un tours de plus !!!! "<<endl;
			for(int j = 0; j< pow(2,i);j++)
			{
				Entiers *Q = new Entiers();
				gggp(g,part[j],Q,part);
				part.push_back(Q);
			}

			/*cout<<"affichage de la partiton en cours : "<<endl;
			for(int k=0; k<part.size(); k++)
			{
				for(int j=0; j<part[k].size(); j++)
				{
					cout<<part[k][j]<<endl;
				}
				cout<<"\n"<<endl;
			}*/
			//cout<<"taille du tableau : "<<part.size()<<endl;
		}
	}
	else
	{
		//cout<<"je jsuis dans gggp_pond"<<endl;

		for(int i = 0; i<floor(log(nbr_parties)/log(2)); i++)
		{
			cout<<"Et un tours de plus !!!! "<<endl;
			for(int j = 0; j< pow(2,i);j++)
			{
				Entiers *Q = new Entiers();
				gggp_pond(g,part.at(j),Q,part);
				clog<<"sortie du gggp_pond"<<endl;
				part.push_back(Q);
			}

			cout<<"****"<<endl;
			for(int k=0; k<part.size(); k++)
			{
				for(int j=0; j<part[k]->size(); j++)
				{
					cout<<part.at(k)->at(j)<<endl;
				}
				cout<<"\n"<<endl;
			}
			cout<<"****"<<endl;
			cout<<"taille du tableau : "<<part.size()<<endl;
		}
	}
}

void bissectionRec(GraphNonOriente *g,EntiersEntiers &Partition,int nbr_parties,const string &nom)
{
	if((nbr_parties&(nbr_parties-1))==0)
	{
		//cout<<"C'est de la forme 2l : "<<nbr_parties<<endl;
		Iter_2l(Partition,nbr_parties,g,nom);
	}
	else
	{
		int puissance_2=0;

		Entiers tailles;

		while(pow(2,puissance_2)<nbr_parties)
			puissance_2++;


		Iter_2l(Partition,pow(2,puissance_2-1),g,nom);

		for(unsigned int i = 0; i< Partition.size() -1 ; i++)
		{
			for(EntiersEntiers::iterator it1 = Partition.begin() + i ; it1!=Partition.end(); it1++)
			{
				if((*it1)->size() > Partition.at(i)->size())
					Partition.at(i)->swap(**it1);
			}
		}
		/*cout<<"****"<<endl;
		for(int k=0; k<Partition.size(); k++)
		{
			for(int j=0; j<Partition[k]->size(); j++)
			{
				cout<<Partition[k]->at(j)<<endl;
			}
			cout<<"\n"<<endl;
		}
		cout<<"****"<<endl;
		cout<<"Et on termine les dernières bissections !!!! "<<endl;*/
		for(int j = 0; j<nbr_parties-pow(2,puissance_2-1);j++)
		{
			Entiers *Q = new Entiers();
			if(nom!="gggp_pond")
				gggp(g,Partition.at(j),Q,Partition);
			else
				gggp_pond(g,Partition.at(j),Q,Partition);
			Partition.push_back(Q);
			//Q.clear();
		}
	}
}

/**
 * Fonction réalisant un partitionnement pseudo aléatoire suivant un voisinage.
 * @param *g : adresse d'un graphe de type boost graphe undirected
 * @param Partition : vecteur contenant des vecteurs d'entiers [tableau contenant les parties de la partition]
 * @param nbr_partie : entier correspondant au nombre de parties voulues pour la partition
 * @return
 */

void Pseudo_random_partitioning(GraphNonOriente *g, EntiersEntiers &Partition, int nbr_parties){
	/*
	 * Principe : distribution des sommets de la première partie en plusieurs autres parties
	 * Le partitionnement étant pseudo aléatoire il n'y a pas de contrainte stricte sur le nombre
	 * de sommets par partie
	 */


	int size = Partition.at(0)->size();
	int cpt_sommets=1;
	int val;
	int cpt;

	if(nbr_parties==size){
		for(uint i = 0; i < nbr_parties;i++){
			if(Partition.at(0)->size()!=1)
					{
						val=rand_fini(0,Partition.at(0)->size()-1);//tirage aléatoire d'un sommets
					}
					else
						val=0;
			int vertex = Partition.at(0)->at(val);
			Entiers *part = new Entiers();
			part->push_back(vertex);// ajout du sommet tiré
			suprim_val(*Partition.at(0),vertex);//suppression du sommet dans la premiere partie


		}
	}
	/*
	 * Boucle sur le nombre de partie à réaliser
	 */
	for(uint i = 0; i < nbr_parties-1;i++){
		if(Partition.at(0)->size()!=1)
		{
			val=rand_fini(0,Partition.at(0)->size()-1);//tirage aléatoire d'un sommets
		}
		else
			val=0;
		int vertex = Partition.at(0)->at(val);
		/*
		 * Initialisation d'un pointeur sur un vecteur d'entier, dans notre cas
		 * la n+1 ième partie de la partition
		 */
		Entiers *part = new Entiers();
		part->push_back(vertex);// ajout du sommet tiré
		suprim_val(*Partition.at(0),vertex);//suppression du sommet dans la premiere partie
		cpt=1;

		/*
		 * Pour chaque element de la nouvelle partie faire
		 */
		for(uint j = 0; j<part->size();j++){
			/*
			 * Détermination des voisins de chacun des sommets de cette nouvelle
			 * partie et ajoue de ce voisin si celui-ci est présent dans la première partie (Partition[0])
			 */
			tie(neighbourIt, neighbourEnd) = adjacent_vertices(part->at(j),*g);
			for (; neighbourIt != neighbourEnd; ++neighbourIt){
				if(In_tab(*Partition.at(0),*neighbourIt)==1){
					cout<<"le voisin déplacé est : "<<*neighbourIt<<endl;
					part->push_back(*neighbourIt);
					cpt_sommets++;
					suprim_val(*Partition.at(0),*neighbourIt);
					cpt++;
				}
				/*
				 * Si le nombre moyen de sommets est atteind dans la partie on sort de la boucle des voisins
				 * Même chose si l'on a rencontré le nombre total de sommets
				 */
				if(cpt==(size/nbr_parties)+1)
					break;
				if(cpt_sommets==size)
					break;
			}

			/*
			 * Même chose
			 */
			if(cpt==(size/nbr_parties)+1)
				break;
			if(cpt_sommets==size)
				break;
		}
		Partition.push_back(part);// ajoue de la nouvelle partie à la partition
		if(cpt_sommets==size)
			break;
	}
}

Graphs Multiniveau(int niveau_contraction, GraphNonOriente *g, GraphOriente *go, int nbr_parties, string type_methode, string choix_affinage, string type_cut,Edges &edge_partie , OutputEdgeList &outputedgeslist, InputEdgeList &inputedgelist, Connections &connections){
	EntiersEntiers Partition;
	Entiers *part = new Entiers();
	Base_Graph baseg;
	baseg.push_back(g);
	ListEntiersEntiers liste_corr;
	uint cpt =0;
	while(num_vertices(*baseg.at(cpt))>niveau_contraction)
	{
		contraction_HEM(baseg.at(cpt),baseg,liste_corr);
		cpt++;
	}

	for(int i =0;i<num_vertices(*baseg.at(baseg.size()-1));i++)
	{
		part->push_back(i);
	}
	Partition.push_back(part);

	bissectionRec(baseg.at(baseg.size()-1),Partition,nbr_parties,type_methode);

	ListEntiersEntiers::iterator lit(liste_corr.end());
	lit--;
	for(uint y =0; y<liste_corr.size();y++){
		projection(Partition,lit);
		double cut = Cut_cluster(Partition,*baseg.at(baseg.size()-2-y),type_cut);
		cout<<"Cout de coupe avant affinage : "<<cut<<endl;

		if(choix_affinage=="charge")
			Affinage_equilibrage_charge(baseg.at(baseg.size()-2-y),Partition);
		else
			Affinage_recherche_locale(baseg.at(baseg.size()-2-y),Partition,cut,type_cut);

		cout<<"Cout de coupe après affinage : "<<cut<<endl;
		lit--;
	}

	Graphs Graphes = Graph_Partition(Partition,go,g,outputedgeslist,inputedgelist,connections);

	for(EntiersEntiers::iterator it = Partition.begin(); it != Partition.end(); it++)
	{
		delete *it;
		*it = NULL;
	}

	for(ListEntiersEntiers::iterator it = liste_corr.begin(); it != liste_corr.end(); it++)
	{
		for(EntiersEntiers::iterator it1 = (*it)->begin(); it1 != (*it)->end(); it1++)
		{
			delete *it1;
			*it1 = NULL;
		}
		delete *it;
		*it = NULL;
	}

	for(Base_Graph::iterator it = baseg.begin(); it != baseg.end(); it++)
	{
		delete *it;
		*it = NULL;
	}

	return Graphes;



}
