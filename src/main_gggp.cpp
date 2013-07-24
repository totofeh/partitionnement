#include "gggp.h"

using namespace std;
using namespace boost;

int main(){
	clock_t t;
	t = clock();  /* Lancement de la mesure */

	srand((unsigned)time(NULL));
	GraphNonOriente *g = new GraphNonOriente();
	GraphOriente *go = new GraphOriente();
	construire_graph(g,go);

	EntiersEntiers Partition;
	Entiers *part = new Entiers();
	Base_Graph baseg;
	baseg.push_back(g);
	ListEntiersEntiers liste_corr;
	uint cpt=0;
	while(num_vertices(*baseg.at(cpt))>4)
	{
		contraction_HEM(baseg.at(cpt),baseg,liste_corr);
		cpt++;
	}

	edge_t e1;
	bool found;
	for(uint i=0;i<baseg.size();i++){
		tie(vertexIt, vertexEnd) = vertices((*baseg.at(i)));
		for (; vertexIt != vertexEnd; ++vertexIt)
		{
			cout << *vertexIt << " est connecté avec ";
			tie(neighbourIt, neighbourEnd) = adjacent_vertices(*vertexIt, (*baseg.at(i)));
			for (; neighbourIt != neighbourEnd; ++neighbourIt){
				cout << *neighbourIt << " ";
				tie(e1,found)=edge(vertex(*vertexIt,*baseg.at(i)),vertex(*neighbourIt,*baseg.at(i)),*baseg.at(i));
				cout << "poids arc : "<<(*baseg.at(i))[e1]._weight<<"\n";
			}
			cout<<" et son poids est de "<< (*baseg.at(i))[*vertexIt]._weight<<endl;
		}
		cout<<"\n"<<endl;
	}


	for(int i =0;i<num_vertices(*baseg.at(baseg.size()-1));i++)
	{
		part->push_back(i);
	}
	Partition.push_back(part);

	bissectionRec(baseg.at(baseg.size()-1),Partition,2,"gggp_pond");
	//Pseudo_random_partitioning(g,Partition,3);
	cout<<"Nombre de parties : "<<Partition.size()<<endl;

	clog<<"Resultat de la partition : "<<endl;
	cout<<"****"<<endl;
	for(uint i = 0; i< Partition.size() ; i++)
	{
		for(uint j = 0 ; j<Partition.at(i)->size() ; j++)
		{
			cout<<(*baseg.at(baseg.size()-1))[Partition.at(i)->at(j)]._index<<endl;
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

	Edges edge_partie;
	OutputEdgeList outputedgeslist(Partition.size());
	InputEdgeList inputedgelist;

	Graphs Graphes = Graph_Partition(Partition,go,g,outputedgeslist,inputedgelist);

	for(uint i = 0; i<Graphes.size();i++){
		tie(vertexIto, vertexEndo) = vertices(*Graphes.at(i));
		for (; vertexIto != vertexEndo; ++vertexIto)
		{
			cout << (*Graphes.at(i))[*vertexIto]._index << " est connecté avec ";
			tie(neighbourIto, neighbourEndo) = adjacent_vertices(*vertexIto,*Graphes.at(i));
			for (; neighbourIto != neighbourEndo; ++neighbourIto)
				cout << (*Graphes.at(i))[*neighbourIto]._index << " ";
			cout<<" et son poids est de "<< (*Graphes.at(i))[*vertexIto]._weight<<endl;
		}
		cout<<"\n"<<endl;
	}

	clog<<"OutputEdgeList : "<<endl;
	for(uint i = 0; i< outputedgeslist.size() ; i++)
	{
		for(uint j = 0; j< outputedgeslist.at(i).size(); j++){
			cout<<outputedgeslist.at(i).at(j).first<<" "<<outputedgeslist.at(i).at(j).second<<endl;
		}
		cout<<"\n"<<endl;
	}

	clog<<"InputEdgeList : "<<endl;
	for(uint i = 0; i< inputedgelist.size() ; i++)
	{
		for(uint j = 0; j< inputedgelist.at(i).size(); j++){
			cout<<inputedgelist.at(i).at(j).first<<" "<<inputedgelist.at(i).at(j).second<<endl;
		}
		cout<<"\n"<<endl;
	}


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

	for(Graphs::iterator it = Graphes.begin(); it != Graphes.end(); it++)
	{
		delete *it;
		*it = NULL;
	}

	delete go;

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
