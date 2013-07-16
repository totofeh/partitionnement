#include "gggp.h"

using namespace std;
using namespace boost;

void gggp(GraphNonOriente &g,Entiers &sommetsSource, Entiers &sommetsDestination,EntiersEntiers &Partition)
{
	int val;
	Entiers sommets_adj;
	if((sommetsSource.size()-1)==0)
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
				gggp(g,*(Partition.at(i)),sommetsDestination,Partition);
			break;
		}
	}
	else
		val=rand_fini(0,sommetsSource.size()-1);//Tirage aléatoire de l'indice du premier sommet entre 0 et taille du tableau -1

	float poids_max=sommetsSource.size()/2.;
	float poids=1;
	Entiers sommets_cut;

	//clog<<"Etape 1 : "<<endl;
	sommetsDestination.push_back(sommetsSource[val]);
	sommetsSource.erase(sommetsSource.begin() + val);

	if(sommetsSource.size()<2)
		return;

	while(poids<poids_max)
	{
//		for(uint i =0; i< sommetsDestination.size();i++){
//			cout<<sommetsDestination.at(i)<<endl;
//		}
		Liste_Voisin(sommetsDestination,sommets_adj,g);
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
				sommets_cut.push_back(Cout_coupe(sommetsDestination,sommets_adj[i],g));
			}
			/*cout<<"\n"<<endl;
			cout<<"Cout de coupe : "<<endl;
			for(int i=0;i<sommets_cut.size();i++)
			{
				cout<<sommets_cut[i]<<endl;
			}*/
			int tmp = recherche_val(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()));
			sommetsDestination.push_back(sommets_adj[tmp]);
			//cout<<"\n"<<endl;
			//cout<<"indice du coup de coupe minimum : "<<recherche_val(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))<<endl;
			//cout<<"Valeur du coup de coupe minimum : "<<sommets_adj[recherche_val(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]<<endl;
			suprim_val(sommetsSource, sommets_adj[tmp]);
			suprim_val(sommets_adj, sommets_adj[tmp]);

			sommets_cut.clear();
			poids++;

		}
	}

	for (int i=0; i<sommetsSource.size();i++)
	{
		for (int j=0; j<sommetsDestination.size();j++)
		{
			remove_edge(sommetsSource[i],sommetsDestination[j],g);
		}
	}

	sort(sommetsDestination);
}

void gggp_pond(GraphNonOriente &g,Entiers &sommetsSource, Entiers &sommetsDestination,EntiersEntiers &Partition){
	property_map<GraphNonOriente,vertex_degree_t>::type poids_sommets=get(vertex_degree_t(),g);
	int val;
	Entiers sommets_adj;
	if((sommetsSource.size()-1)==0){
		val=0;
		cout<<"Entré dans le debug ! "<<endl;
		Entiers tailles;
		for(int i=0;i<Partition.size();i++){
			tailles.push_back((*Partition.at(i)).size());
		}
		int tmp=*max_element(tailles.begin(),tailles.end());
		for(int i=0; i<Partition.size();i++){
			if(Partition.at(i)->size()==tmp)
				gggp_pond(g,*Partition[i],sommetsDestination,Partition);
			break;
		}
	}
	else
		val=rand_fini(0,sommetsSource.size()-1);//Tirage aléatoire de l'indice du premier sommet entre 0 et taille du tableau -1
	//cout<<"Sommet de départ : "<<sommetsSource[val]<<endl;
	//	cout<<"\n"<<endl;
	float poids_max=0;
	for(int i=0;i<sommetsSource.size();i++){
		poids_max+=poids_sommets[sommetsSource[i]];
	}
	poids_max/=2.;
	cout<<"Le moitié du poids total est de : "<<poids_max<<endl;
	cout<<"\n"<<endl;
	float poids=poids_sommets[sommetsSource[val]];
	vector<float> sommets_cut;

	sommetsDestination.push_back(sommetsSource.at(val));
	sommetsSource.erase(sommetsSource.begin() + val);

	if(sommetsSource.size()<2)
			return;

	while(poids<poids_max)
	{
		Liste_Voisin(sommetsDestination,sommets_adj,g);
		if((sommets_adj.size()==0))
		{
			cout<<"Je suis sorti !!!! "<<endl;
			return;
		}
		else{
			sort(sommets_adj);

			for(int i=0;i<sommets_adj.size();i++)
			{
				sommets_cut.push_back(Cout_coupe_pond(sommetsDestination,sommets_adj[i],g));
			}
			sommetsDestination.push_back(sommets_adj[recherche_val2(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]);
			poids+=poids_sommets[sommets_adj[recherche_val2(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]];
			suprim_val(sommetsSource, sommets_adj[recherche_val2(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]);
			suprim_val(sommets_adj, sommets_adj[recherche_val2(sommets_cut,*min_element(sommets_cut.begin(),sommets_cut.end()))]);

			sommets_cut.clear();

		}
	}

	for (int i=0; i<sommetsSource.size();i++)
	{
		for (int j=0; j<sommetsDestination.size();j++)
		{
			remove_edge(sommetsSource[i],sommetsDestination[j],g);
		}
	}

	sort(sommetsDestination);
}

void Iter_2l(EntiersEntiers &part, int nbr_parties, GraphNonOriente &g,const string &nom)
{
	if (nom=="gggp"){
		//cout<<"je jsuis dans gggp"<<endl;

		for(int i = 0; i<floor(log(nbr_parties)/log(2)); i++)
		{
			//cout<<"Et un tours de plus !!!! "<<endl;
			for(int j = 0; j< pow(2,i);j++)
			{
				Entiers *Q = new Entiers();
				gggp(g,*part[j],*Q,part);
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
				gggp_pond(g,*part[j],*Q,part);
				part.push_back(Q);
			}

			/*cout<<"****"<<endl;
			for(int k=0; k<part.size(); k++)
			{
				for(int j=0; j<part[k]->size(); j++)
				{
					cout<<part.at(k)->at(j)<<endl;
				}
				cout<<"\n"<<endl;
			}
			cout<<"****"<<endl;
			cout<<"taille du tableau : "<<part.size()<<endl;*/
		}
	}
}

void bissectionRec(GraphNonOriente &g,EntiersEntiers &Partition,int nbr_parties,const string &nom)
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
			Entiers * Q = new Entiers();
			if(nom=="gggp")
				gggp(g,*Partition.at(j),*Q,Partition);
			else
				gggp_pond(g,*Partition.at(j),*Q,Partition);
			Partition.push_back(Q);
			//Q.clear();
		}
	}
}
