#include "gggp.h"

using namespace std;
using namespace boost;

int main(){
	clock_t t;
	t = clock();  /* Lancement de la mesure */

	srand((unsigned)time(NULL));
	GraphNonOriente *g = new GraphNonOriente();
	construire_graph(g);

	EntiersEntiers Partition;
	Entiers *part = new Entiers();
	Base_Graph baseg;
	baseg.push_back(g);
	ListEntiersEntiers liste_corr;
	uint cpt=0;
	//while(num_vertices(*baseg.at(cpt))>4)
	//{
		contraction_HEM(baseg.at(cpt),baseg,liste_corr);
		cpt++;
	//}

	edge_t e1;
	bool found;
	for(uint i=0;i<baseg.size();i++){
		property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),*baseg.at(i));
		tie(vertexIt, vertexEnd) = vertices((*baseg.at(i)));
		property_map<GraphNonOriente,vertex_degree_t>::type poids_sommets=get(vertex_degree_t(),(*baseg.at(i)));
		for (; vertexIt != vertexEnd; ++vertexIt)
		{
			cout << *vertexIt << " est connecté avec ";
			tie(neighbourIt, neighbourEnd) = adjacent_vertices(*vertexIt, (*baseg.at(i)));
			for (; neighbourIt != neighbourEnd; ++neighbourIt){
				cout << *neighbourIt << " ";
				tie(e1,found)=edge(vertex(*vertexIt,*baseg.at(i)),vertex(*neighbourIt,*baseg.at(i)),*baseg.at(i));
				cout << "poids arc : "<<get(poids_arc,e1)<<"\n";
			}
			cout<<" et son poids est de "<< poids_sommets[*vertexIt]<<endl;
		}
		cout<<"\n"<<endl;
	}


	for(int i =0;i<num_vertices(*(baseg.at(baseg.size()-1)));i++)
	{
		part->push_back(i);
	}
	Partition.push_back(part);

	bissectionRec(*(baseg.at(baseg.size()-1)),Partition,3,"gggp_pond");
	cout<<"Nombre de parties : "<<Partition.size()<<endl;

	clog<<"Resultat de la partition : "<<endl;
	cout<<"****"<<endl;
	for(uint i = 0; i< Partition.size() ; i++)
	{
		for(uint j = 0 ; j<Partition.at(i)->size() ; j++)
		{
			cout<<Partition.at(i)->at(j)<<endl;
		}
		cout<<"\n"<<endl;
	}
	cout<<"****"<<endl;


	ListEntiersEntiers::iterator lit(liste_corr.end());
	lit--;
	for(uint y =0; y<liste_corr.size();y++){
		projection(Partition,lit);

		clog<<"liste de correspondance : "<<endl;
		for(uint i = 0; i < (*lit)->size(); i++)
		{
			for(uint j = 0; j < (*lit)->at(i)->size();j++){
				cout<<(*lit)->at(i)->at(j)<<endl;;
			}
			cout<<"\n"<<endl;
		}

		clog<<"Resultat projection : "<<endl;
		for(uint i = 0; i< Partition.size() ; i++)
			{
				for(uint j = 0 ; j<Partition.at(i)->size() ; j++)
				{
					cout<<Partition.at(i)->at(j)<<endl;
				}
				cout<<"\n"<<endl;
			}

		double cut = Cut_cluster(Partition,*baseg.at(baseg.size()-2-y),"norm");

		cout<<"Cout de coupe avant affinage : "<<cut<<endl;

		Affinage_recherche_locale(baseg.at(baseg.size()-2-y),Partition,cut,"norm");
		//Affinage_equilibrage_charge(baseg.at(baseg.size()-2-y),Partition);
		cout<<"Cout de coupe après affinage : "<<cut<<endl;
		cout<<"\n"<<endl;
		double tmp = Cut_cluster(Partition,*baseg.at(baseg.size()-2-y),"norm");
		cout<<"verification cout de coupe après affinage : "<<tmp<<endl;
		cout<<"\n"<<endl;
		clog<<"Partition après affinage : "<<endl;
		for(uint i = 0; i< Partition.size() ; i++)
		{
			for(uint j = 0 ; j<Partition.at(i)->size() ; j++)
			{
				cout<<Partition.at(i)->at(j)<<endl;
			}
			cout<<"\n"<<endl;
		}

		lit--;
	}
	cout<<"mathieu va me buter ! et en plus c'est walker !!!"<<endl;
	cout<<"\n"<<endl;

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


	 t = clock() - t;
	 printf ("It took me (%f seconds).\n",((double)t)/CLOCKS_PER_SEC);

	 //EntiersEntiersEntiers Stock_Partition;
	/*for(int i=0;i<11;i++){
		g1=g;
		Partition1=Partition;
		bissectionRec(g1,Partition1,4,"gggp");
		Stock_Partition.push_back(Partition1);
		Cut.push_back(Cut_cluster(Partition1,g));
		Partition1.clear();
		g1.clear();
	}

	for(int i =0;i<Cut.size();i++){
		cout<<Cut[i]<<endl;
	}
	cout<<"\n"<<endl;
	cout<<recherche_val_double(Cut,*min_element(Cut.begin(),Cut.end()))<<endl;
	cout<<"\n"<<endl;
	EntiersEntiers tmp = Stock_Partition[recherche_val_double(Cut,*min_element(Cut.begin(),Cut.end()))];
	for(int i =0; i<tmp.size();i++){
		for(int j=0; j<tmp[i].size();j++){
			cout<<tmp[i][j]<<endl;
		}
		cout<<"\n"<<endl;
	}*/

	/*EntiersEntiers liste_corr;
	//double cut = Cut_cluster(Partition,g);
	Base_Graph baseg;
	baseg.push_back(g);




	cout<<"LIste des noeuds fusionés : "<<endl;
	for(uint i = 0; i< liste_corr.size() ; i++)
		{
			for(uint j = 0 ; j<liste_corr.at(i).size() ; j++)
			{
				cout<<liste_corr[i][j]<<endl;
			}
			cout<<"\n"<<endl;
		}*/


};
