#include "util.h"

using namespace std;
using namespace boost;

GraphNonOriente::vertex_iterator vertexIt, vertexEnd;
GraphNonOriente::adjacency_iterator neighbourIt, neighbourEnd;

struct myclass{
	bool operator() (Entiers *i , Entiers *j) {return (i->at(0)<j->at(0));}
}myobject;

struct myclass2{
	bool operator() (Entiers *i , Entiers *j, GraphNonOriente *g) {return (Calcul_poids(i,g)<Calcul_poids(j,g));}
}mon_tri;

bool Est_connexe(GraphNonOriente *g, EntiersEntiers Partition, Entiers &part){
	GraphNonOriente copie_g;
	copie_g = *g;

	for (uint i=0; i<Partition.size()-1;i++)
	{
		for (uint j=1+i; j<Partition.size();j++)
		{
			for (uint k=0; k<Partition.at(i)->size();k++)
			{
				for (uint y=0; y<Partition.at(j)->size();y++)
				{
					remove_edge(Partition.at(i)->at(k),Partition.at(j)->at(y),copie_g);
				}
			}
		}
	}

	int val;
	Entiers sommets;

	if(part.size()==1)
		val = 0;
	else
		val=rand_fini(0,part.size()-1);

	int vertex = part.at(val);
	sommets.push_back(vertex);

	for(uint i = 0;i<sommets.size();i++){
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(sommets.at(i),copie_g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt){
			if(In_tab(sommets,*neighbourIt)!=1)
				sommets.push_back(*neighbourIt);
		}
	}

	if(part.size()!=sommets.size())
		return false;
	else
		return true;

}

void projection(EntiersEntiers &Partition,ListEntiersEntiers::iterator lit){

		EntiersEntiers new_partition;
		for(uint i = 0; i< Partition.size() ; i++)
		{
			Entiers *new_part = new Entiers();
			for(uint j = 0 ; j<Partition.at(i)->size() ; j++)
			{
				for(uint k = 0; k<((*lit)->at(Partition.at(i)->at(j)))->size();k++){
					new_part->push_back((*lit)->at(Partition.at(i)->at(j))->at(k));
				}

			}
			new_partition.push_back(new_part);
		}

		for(EntiersEntiers::iterator it = Partition.begin(); it != Partition.end(); it++)
		{
			delete *it;
			*it = NULL;
		}

		Partition = new_partition;
		for(uint i =0; i<Partition.size(); i++){
			sort(*Partition.at(i));
		}

		new_partition.clear();
}

double Calcul_poids(Entiers *partie, GraphNonOriente *g){
	property_map<GraphNonOriente,vertex_degree_t>::type poids_sommets=get(vertex_degree_t(),*g);
	double poids=0;

	for(uint j = 0; j<partie->size();j++){
		poids+=poids_sommets[partie->at(j)];
	}

	return poids;
}

void Affinage_equilibrage_charge(GraphNonOriente *g, EntiersEntiers &Partition){
	property_map<GraphNonOriente,vertex_degree_t>::type poids_sommets=get(vertex_degree_t(),*g);
	/*
	 * Calcule de la somme des poids du graphe et le poids moyen à atteindre
	 */
	double poids_moyen = 0.;
	for(uint i =0; i<num_vertices(*g); i++)
		poids_moyen+=poids_sommets[i];

	poids_moyen /=Partition.size();
	vector<double> poids_parties;

	Entiers valeurs_interdites;

	/*
	 * Calcule du poids de chaque partie de la partition
	 */
	for(uint i = 0; i<Partition.size();i++){
		double tmp = Calcul_poids(Partition.at(i),g);
		poids_parties.push_back(tmp);
	}

	clog<<"Poids initial des parties : "<<endl;
	for(uint i = 0; i<poids_parties.size(); i++){
		cout<<poids_parties.at(i)<<" ";
	}
	cout<<"\n"<<endl;

	double critere = 0.;

	for(uint i = 0; i<poids_parties.size();i++){
		critere += abs(poids_parties.at(i)-poids_moyen);
	}
	critere/=Partition.size();

	double p_max = *max_element(poids_parties.begin(),poids_parties.end());

	cout<<"Valeurs du criètre de départ : "<<critere<<endl;
	double best_critere = critere-0.0000001;
	int nbr_passage = 1;

	while(best_critere<critere || nbr_passage<Partition.size()){

		critere = best_critere;
		int cpt = recherche_val_double(poids_parties,p_max);
		cout<<"Partie de poids maximum : "<<cpt<<endl;

		bool decision = false;
		int nbr_pass_interne = 0;

		Entiers random_orders(Partition.at(cpt)->size());
		for (uint i=0 ; i<Partition.at(cpt)->size() ; i++)
			random_orders.at(i)=Partition.at(cpt)->at(i);

		for (uint j=0 ; j<Partition.at(cpt)->size()-1 ; j++) {
			int rand_pos = (rand() % Partition.at(cpt)->size()-j)+j;
			int tmp = random_orders[j];
			random_orders[j] = random_orders[rand_pos];
			random_orders[rand_pos] = tmp;
		}

		int size;

		if(Partition.at(cpt)->size()>400)
			size = 400;
		else
			size = Partition.at(cpt)->size();

		while(decision!=true && nbr_pass_interne < size){
			int vertex = random_orders.at(nbr_pass_interne);
			Entiers community = Neigh_community(g,Partition,vertex,cpt);
			if(community.size()!=0)
			{
				vector<double> new_poids_parties = poids_parties;
				vector<double> tmp_critere;
				for(uint k = 0; k < community.size();k++)
				{
					new_poids_parties = poids_parties;
					new_poids_parties.at(community.at(k))+=poids_sommets[vertex];
					new_poids_parties.at(cpt)-=poids_sommets[vertex];

					double new_critere = 0.;

					for(uint i = 0; i<poids_parties.size();i++){
						new_critere += abs(new_poids_parties.at(i)-poids_moyen);
					}
					new_critere/=Partition.size();
					tmp_critere.push_back(new_critere);
				}
				double val_min = *min_element(tmp_critere.begin(),tmp_critere.end());
				int val = recherche_val_double(tmp_critere,val_min);
				if(val_min<critere && poids_parties.at(val)!=-1)
				{
					Partition.at(community.at(val))->push_back(vertex);
					suprim_val(*Partition.at(cpt),vertex);
					sort(*Partition.at(community.at(val)));

					if(Est_connexe(g,Partition,*Partition.at(cpt))!=1)
					{
						suprim_val(*Partition.at(community.at(val)),vertex);
						Partition.at(cpt)->push_back(vertex);
						sort(*Partition.at(cpt));
						cout<<" C'EST MORT RETOUR EN ARRIERE ! "<<endl;
					}
					else
					{
						poids_parties = new_poids_parties;
						decision = true;
						cout<<" Modification reussi ! "<<endl;
					}
				}
			}
			nbr_pass_interne++;
		}
		if(decision==false)
		{
			nbr_passage++;
			valeurs_interdites.push_back(cpt);
			poids_parties.at(cpt)=-1;
			clog<<"nouveau passag ! "<<endl;
		}
		else
		{
			best_critere = 0.;

			for(uint i = 0; i<poids_parties.size();i++){
				best_critere += abs(poids_parties.at(i)-poids_moyen);
			}
			best_critere/=Partition.size();
			nbr_passage = 0;
		}

		clog<<"Poids des parties modifié : "<<endl;
		for(uint i = 0; i<poids_parties.size(); i++){
			cout<<poids_parties.at(i)<<" ";
		}
		cout<<"\n"<<endl;
		p_max = *max_element(poids_parties.begin(),poids_parties.end());
		cout<<"Valeurs du criètre : "<<critere<<endl;
		cout<<"Valeurs du best_criètre : "<<best_critere<<endl;
		cout<<"Nombre de passage : "<<nbr_passage<<endl;
		cout<<"\n"<<endl;

	}
}

void Affinage_recherche_locale(GraphNonOriente *g, EntiersEntiers &Partition, double &cut, string name){

	Entiers random_orders(num_vertices(*g)); //gestion d'un tableau contenant tout les sommets et ranger de façon aléatoire

	for (uint i=0 ; i<random_orders.size() ; i++)
		random_orders.at(i)=i;

	for (uint j=0 ; j<num_vertices(*g)-1 ; j++) {
		int rand_pos = (rand() % num_vertices(*g)-j)+j;
		int tmp = random_orders[j];
		random_orders[j] = random_orders[rand_pos];
		random_orders[rand_pos] = tmp;
	}
	int size = random_orders.size();

	if(random_orders.size()>400)
		size=400;

	for(uint i = 0; i<size;i++){
		if(random_orders.at(i)!=-1){
			int vertex = random_orders.at(i);

			int comm = In_community_dichotomie(Partition, vertex);
			Entiers community = Neigh_community(g,Partition,vertex,comm);
			vector<double> tmp_cut;

			if(community.size()!=0 && Partition.at(comm)->size()!=1){
				tmp_cut = modif_cut_tmp(g,Partition,vertex,comm,community,cut,name);
				for(uint z = 0; z<tmp_cut.size(); z++){
					cout<<tmp_cut.at(z)<<endl;
				}
				cout<<"\n"<<endl;
				double cut_min = *min_element(tmp_cut.begin(),tmp_cut.end());
				cout<<"cout de coupe minimum de la liste : "<<cut_min<<endl;
				if(cut_min<cut){
					clog<<"Changement ! "<<endl;
					int tmp = recherche_val_double(tmp_cut,cut_min);
					cut=cut_min;
					Partition.at(community.at(tmp))->push_back(vertex);
					suprim_val(*Partition.at(comm),vertex);
					sort(*Partition.at(community.at(tmp)));

				}
			}

				//Modif_fonction_Gain_Cut(Partition,g,community,vertex,cut,name);


			/*if(Est_connexe(g,Partition,*Partition.at(comm))!=1)
			{
				suprim_val(*Partition.at(community.at(tmp)),vertex);
				Partition.at(comm)->push_back(vertex);
				sort(*Partition.at(comm));
				cout<<" C'EST MORT RETOUR EN ARRIERE ! "<<endl;
			}*/

		}
	}
}

Entiers Neigh_community(GraphNonOriente *g, EntiersEntiers &Partition, int vertex, int comm_in){
	Entiers community;
	int comm;
	tie(neighbourIt, neighbourEnd) = adjacent_vertices(vertex,*g);
	for (; neighbourIt != neighbourEnd; ++neighbourIt){
		comm = In_community_dichotomie(Partition,*neighbourIt);
		if(In_tab(community,comm)!=1 && comm!=comm_in)
			community.push_back(comm);
	}
	return community;
}

double Modif_Cut_one_cluster(Entiers &cluster, GraphNonOriente &g, double &vol)
{
	property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),g);
	edge_t e1;
	bool found;
	double cpt=0.;
	for(int i=0;i<cluster.size();i++){
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(cluster[i], g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt){
			tie(e1,found)=edge(vertex(cluster[i],g),vertex(*neighbourIt,g),g);
			if(In_tab(cluster,*neighbourIt)!=1){
				cpt+=get(poids_arc,e1);
			}
		}
	}
	vol = Cluster_Degree(g,cluster);
	return cpt/2.;

}

vector<double> modif_cut_tmp(GraphNonOriente *g, EntiersEntiers &Partition, int vertexs, int comm_in, Entiers community, double cut,string name){
	property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),*g);
	edge_t e1;
	bool found;
	cout<<"le sommet tiré est : "<<vertexs<<endl;

	if(name!="norm")
	{
		vector<double> modif_cut(community.size());
		double cpt_comm_in;
		double cpt_comm_out;
		for(uint i =0; i<community.size(); i++){
			double tmp_cut = cut;
			cpt_comm_in=0;
			cpt_comm_out=0;

			tie(neighbourIt, neighbourEnd) = adjacent_vertices(vertexs,*g);
			for (; neighbourIt != neighbourEnd; ++neighbourIt){
				tie(e1,found)=edge(vertex(vertexs,*g),vertex(*neighbourIt,*g),*g);
				if(In_tab(*Partition.at(comm_in),*neighbourIt)==1)
					cpt_comm_in+=get(poids_arc,e1);
				else if(In_tab(*Partition.at(community.at(i)),*neighbourIt)==1)
					cpt_comm_out+=get(poids_arc,e1);
			}

			tmp_cut+=cpt_comm_in;
			tmp_cut-=cpt_comm_out;

			modif_cut.at(i)=tmp_cut;
		}
		return modif_cut;
	}
	else{
		vector<double> modif_cut(community.size());
		double cpt_comm_in;
		double cpt_comm_out;
		double tmp_cut;

		for(uint i =0; i<community.size(); i++){
			vector<vector<double> > tab_cut;
			for(uint k = 0; k < Partition.size();k++){
				vector<double> tmp;
				double vol = 0.;
				double cut = Modif_Cut_one_cluster(*Partition.at(k), *g, vol);
				tmp.push_back(cut);
				tmp.push_back(vol);
				tab_cut.push_back(tmp);
			}

			tmp_cut =0.;
			cpt_comm_in=0.;
			cpt_comm_out=0.;

			tie(neighbourIt, neighbourEnd) = adjacent_vertices(vertexs,*g);
			for (; neighbourIt != neighbourEnd; ++neighbourIt){
				tie(e1,found)=edge(vertex(vertexs,*g),vertex(*neighbourIt,*g),*g);
				if(In_tab(*Partition.at(comm_in),*neighbourIt)==1)
					cpt_comm_in+=get(poids_arc,e1);
				else if(In_tab(*Partition.at(community.at(i)),*neighbourIt)==1)
					cpt_comm_out+=get(poids_arc,e1);
			}

			cpt_comm_in/=2.;
			cpt_comm_out/=2.;

			tab_cut.at(comm_in).at(0)+=cpt_comm_in;
			tab_cut.at(comm_in).at(0)-=cpt_comm_out;
			tab_cut.at(comm_in).at(1)-= Degree(*g ,vertexs);

			tab_cut.at(community.at(i)).at(0)+=cpt_comm_in;
			tab_cut.at(community.at(i)).at(0)-=cpt_comm_out;
			tab_cut.at(community.at(i)).at(1)+= Degree(*g ,vertexs);

			for(uint j = 0; j < tab_cut.size();j++){
				tmp_cut+=((tab_cut.at(j).at(0))/(tab_cut.at(j).at(1)));
			}

			modif_cut.at(i)=tmp_cut;
		}
		return modif_cut;
	}


}

void Modif_fonction_Gain_Cut(EntiersEntiers &Partition,GraphNonOriente *g, Entiers &community, int val, double &cut,string name)
{
	/*cout<<"Nombre de communauté voisine : "<<community.size()<<endl;
	cout<<"\n"<<endl;*/
	for(uint i = 0; i<community.size();i++){
		EntiersEntiers new_partition;
		for(uint k = 0; k < Partition.size();k++){
			Entiers * tmp = new Entiers();
			for(uint j = 0;j<Partition.at(k)->size();j++){
				tmp->push_back(Partition.at(k)->at(j));
			}
			new_partition.push_back(tmp);
		}

		/*cout<<"Avant Modification partition"<<endl;
		cout<<"************"<<endl;
		for(uint t = 0; t< new_partition.size() ; t++)
				{
					for(uint j = 0 ; j<new_partition.at(t)->size() ; j++)
					{
						cout<<new_partition.at(t)->at(j)<<endl;
					}
					cout<<"\n"<<endl;
				}
		cout<<"************"<<endl;*/


		new_partition.at(community.at(i))->push_back(val);
		suprim_val(*new_partition.at(In_community_dichotomie(Partition,val)),val);
		sort(*new_partition.at(community.at(i)));

		/*cout<<"Modification partition"<<endl;
		cout<<"************"<<endl;
		for(uint t= 0; t< new_partition.size() ; t++)
		{
			for(uint j = 0 ; j<new_partition.at(t)->size() ; j++)
			{
				cout<<new_partition.at(t)->at(j)<<endl;
			}
			cout<<"\n"<<endl;
		}
		cout<<"************"<<endl;*/

		double coupe = Cut_cluster(new_partition,*g,name);

		//cout<<"cout de coupe : "<<coupe<<endl;
		if(coupe<cut)
		{
			for(EntiersEntiers::iterator it = Partition.begin(); it != Partition.end(); it++)
			{
				delete *it;
				*it = NULL;
			}
			Partition = new_partition;
			cut = coupe;
		}
		else
		{
			for(EntiersEntiers::iterator it = new_partition.begin(); it != new_partition.end(); it++)
			{
				delete *it;
				*it = NULL;
			}
		}
	}
}

void contraction_HEM(GraphNonOriente *g, Base_Graph &baseg, ListEntiersEntiers &liste_corr){
	GraphNonOriente *gtmp = new GraphNonOriente();
	*gtmp=*g;
	property_map<GraphNonOriente,vertex_degree_t>::type poids_sommets=get(vertex_degree_t(),*gtmp);
	property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),*gtmp);
	Entiers Random_list_vertices; // Initialisation du tableau de sommets rangés aléatoirements
	EntiersEntiers *tableau_de_correspondance = new EntiersEntiers();
	edge_t e1,e2; // Iterateurs sur les arcs
	bool found;
	int nbr_vertex = num_vertices(*gtmp);
	Entiers sommets_a_detruire; // Initialisation d'un tableau pret à recevoir les "sommets à détruire"
	/*
	 * Création d'un vecteur contenant l'ensemble des sommets du graphe. Ces sommets sont rangés
	 * aléatoirement afin de simuler un tirage aléatoire
	 */
	for (uint i=0 ; i<nbr_vertex ; i++)
		Random_list_vertices.push_back(i);
	for (uint j=0 ; j<nbr_vertex-1 ; j++) {
		int rand_pos = rand()%(nbr_vertex-j)+j;
		int tmp      = Random_list_vertices[j];
		Random_list_vertices[j] = Random_list_vertices[rand_pos];
		Random_list_vertices[rand_pos] = tmp;
	}

	/*
	 * Pour chaque sommet non verrouiller faire ....
	 */
	for(uint i=0; i<nbr_vertex; i++){
		int vertexs = Random_list_vertices[i];
		if(vertexs!=-1){
			Entiers liste_voisin = Liste_adjacence(*gtmp,vertexs,Random_list_vertices); // Recherche des sommets adjacents au sommets  tiré
			if(liste_voisin.size()!=0){
				/*
				 * S'il en existe au mois un sommet adjacent au sommet tiré qui n'est pas verrouillé, on
				 * choisi celui dont l'arc les reliants est le plus fort. Dans le cas où les arcs ont tous
				 * le même poids, on selectionne le sommet d'identifiant le plus petit
				 */
				double poids_a=0.;
				int best_vertexs;
				for(uint j=0;j<liste_voisin.size();j++){
					tie(e1,found)=edge(vertex(vertexs,*gtmp),vertex(liste_voisin[j],*gtmp),*gtmp);
					if(get(poids_arc,e1)>poids_a){
						best_vertexs = liste_voisin[j];
						poids_a = get(poids_arc,e1);
					}
				}


				Entiers * couple = new Entiers(); // Initialisation du vecteur contenant le couple de sommet fusionné
				int vertex_delete = max(vertexs,best_vertexs); // Sommet d'indentifiant le plus grand (qui sera détruit)
				//cout<<"sommet détruit : "<<vertex_delete<<endl;
				int vertex_save = min(vertexs,best_vertexs); // Sommet d'identifiant le plus petit (qui sera conservé)
				//cout<<"sommet sauvé : "<<vertex_save<<endl;

				sommets_a_detruire.push_back(vertex_delete); // On ajoute le sommet détruit au tableau des sommets à détruire
				/*
				 * On ajoute au tableau "couple" le couple de sommet à fusionner
				 */
				couple->push_back(vertex_save);
				couple->push_back(vertex_delete);
				tableau_de_correspondance->push_back(couple); // Ajout du "couple" à la liste de correspondance

				remove_edge(vertex_save,vertex_delete,*gtmp); // Suppression de l'arc reliant le couple de sommets

				Entiers neigh_vertex_save; // Initialisation du vecteur contenant les somemts adjacents au "sommet sauvegardé"
				Entiers neigh_vertex_delete; // Initialisation du vecteur contenant les somemts adjacents au "sommet à détruire"
				tie(neighbourIt, neighbourEnd) = adjacent_vertices(vertex_save,*gtmp);

				/*
				 * Remplissage de ces deux tableaux à l'aide de la fonction adjacent_vertices de boost graph
				 * [La création de ces tableaux est nécéssaire du fait que certains arcs sont détruit au cours
				 * du processus]
				 */
				for (; neighbourIt != neighbourEnd; ++neighbourIt){
					neigh_vertex_save.push_back(*neighbourIt);
				}

				tie(neighbourIt, neighbourEnd) = adjacent_vertices(vertex_delete,*gtmp);
				for (; neighbourIt != neighbourEnd; ++neighbourIt){
					neigh_vertex_delete.push_back(*neighbourIt);
				}

				/*
				 * Recherche de sommets communs entre le "sommet sauvegardé" et le "sommet à détruire"
				 * S'il existe un tel sommet "v" alors on ajoute le poids de l'arcs (vertex_delet,v)
				 * à celui de l'arcs (vertex_save,v) et on détruit l'arcs reliant "v" au "sommet à détruire"
				 */
				for(int j=0;j<neigh_vertex_delete.size();j++){
					if(In_tab(neigh_vertex_save,neigh_vertex_delete[j])==1){
						tie(e2,found)=edge(vertex(vertex_save,*gtmp),vertex(neigh_vertex_delete[j],*gtmp),*gtmp);
						tie(e1,found)=edge(vertex(vertex_delete,*gtmp),vertex(neigh_vertex_delete[j],*gtmp),*gtmp);
						get(poids_arc,e2)+=get(poids_arc,e1);
						remove_edge(vertex_delete,neigh_vertex_delete[j],*gtmp);
					}
					else
					{
						tie(e1,found)=edge(vertex(vertex_delete,*gtmp),vertex(neigh_vertex_delete[j],*gtmp),*gtmp);
						add_edge(vertex_save,neigh_vertex_delete[j],get(poids_arc,e1),*gtmp);
						remove_edge(vertex_delete,neigh_vertex_delete[j],*gtmp);
					}
				}

				poids_sommets[vertex_save]+=poids_sommets[vertex_delete]; // ajout du poids du sommet détruit au sommet conservé
				/*
				 * Vérouillage du "sommet sauvegardé" et du "sommet à détruire"
				 */
				Random_list_vertices[i]=-1;
				Random_list_vertices[recherche_val(Random_list_vertices,best_vertexs)]=-1;
			}
			else{
				/*
				 * Et si le sommet tiré ne possède pas de sommet adjacent non verrouillé
				 * alors on l'ajoute à la liste de correspondance des sommets et on
				 * le verrouille
				 */
				Entiers *couple = new Entiers();
				couple->push_back(Random_list_vertices.at(i));
				tableau_de_correspondance->push_back(couple);
				Random_list_vertices[i]=-1;
			}
		}
	}

	sort(sommets_a_detruire); // Trie dans l'ordre croissant des "sommets à détruire"
	//cout<<"\n"<<endl;
	/*
	 * Suppression des sommets de la liste "sommets à détruire". Cette suppression est
	 * effectuée dans l'ordre décroissant afin à maintenir à jour la renumérotation
	 * des sommets
	 */
	for(int j=(sommets_a_detruire.size()-1);j>-1;j--){
		//cout<<"Noeuds a supprimer : "<<sommets_a_detruire.at(j)<<endl;
		remove_vertex(sommets_a_detruire[j],*gtmp);
	}

	/**clog<<"Affichage avant tri "<<endl;
	for(uint k = 0;k<tableau_de_correspondance->size();k++){
		for(uint v = 0; v<tableau_de_correspondance->at(k)->size();v++){
			cout<<tableau_de_correspondance->at(k)->at(v)<<" ";
		}
		cout<<"\n"<<endl;
	}*/
	sort(tableau_de_correspondance->begin(),tableau_de_correspondance->end(),myobject); // Trie dans l'ordre croissant des couples de sommets de la liste de correspondance

	clog<<"Tableau de correspondance "<<endl;
	for(uint k = 0;k<tableau_de_correspondance->size();k++){
		for(uint v = 0; v<tableau_de_correspondance->at(k)->size();v++){
			cout<<tableau_de_correspondance->at(k)->at(v)<<" ";
		}
		cout<<"\n"<<endl;
	}

	liste_corr.push_back(tableau_de_correspondance);
	cout<<"\n"<<endl;
	baseg.push_back(gtmp); // Ajout du graphe modifié à la "base des graphe"
}

Entiers Liste_adjacence(GraphNonOriente &g, int vertexs,const Entiers &random_vertices){ // a revoir !!!!
	Entiers liste_voisin;
	tie(neighbourIt, neighbourEnd) = adjacent_vertices(vertexs, g);
	for (; neighbourIt != neighbourEnd; ++neighbourIt){
		if(In_tab(random_vertices,*neighbourIt)==1)
			liste_voisin.push_back(*neighbourIt);
	}
	return liste_voisin;
}

void construire_graph(GraphNonOriente *g){
	property_map<GraphNonOriente,vertex_degree_t>::type poids_sommets=get(vertex_degree_t(),*g);
	add_edge(0,1,1,*g);
	add_edge(0,2,1,*g);
	add_edge(0,3,1,*g);
	add_edge(1,2,1,*g);
	add_edge(1,4,1,*g);
	add_edge(1,5,1,*g);
	add_edge(1,6,1,*g);
	add_edge(2,6,1,*g);
	add_edge(2,3,1,*g);
	add_edge(3,9,1,*g);
	add_edge(3,10,1,*g);
	add_edge(4,5,1,*g);
	add_edge(5,6,1,*g);
	add_edge(4,7,1,*g);
	add_edge(4,8,1,*g);
	add_edge(7,8,1,*g);
	add_edge(9,10,1,*g);

	put(poids_sommets,0,3);

	for(int i=1;i<4;i++)
		put(poids_sommets,i,2);
	for(int i=4;i<7;i++)
		put(poids_sommets,i,1.5);
	put(poids_sommets,7,1);
	put(poids_sommets,8,1);
	put(poids_sommets,9,1.5);
	put(poids_sommets,10,1.5);

	tie(vertexIt, vertexEnd) = vertices(*g);
	for (; vertexIt != vertexEnd; ++vertexIt)
	{
		cout << *vertexIt << " est connecté avec ";
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(*vertexIt, *g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt)
			cout << *neighbourIt << " ";
		cout<<" et son poids est de "<< poids_sommets[*vertexIt];
		cout << "\n";
	}

	/*property_map<GraphNonOriente,vertex_degree_t>::type poids_sommets=get(vertex_degree_t(),*g);
		add_edge(0,1,1,*g);
		add_edge(0,2,1,*g);
		add_edge(0,3,1,*g);
		add_edge(1,2,1,*g);
		add_edge(1,4,1,*g);
		add_edge(1,5,1,*g);
		add_edge(1,6,1,*g);
		add_edge(2,6,1,*g);
		add_edge(2,3,1,*g);
		add_edge(3,18,1,*g);
		add_edge(3,21,1,*g);
		add_edge(4,5,1,*g);
		add_edge(4,7,1,*g);
		add_edge(4,8,1,*g);
		add_edge(4,9,1,*g);
		add_edge(5,6,1,*g);
		add_edge(6,15,1,*g);
		add_edge(6,16,1,*g);
		add_edge(6,17,1,*g);
		add_edge(7,8,1,*g);
		add_edge(7,10,1,*g);
		add_edge(7,11,1,*g);
		add_edge(8,9,1,*g);
		add_edge(8,11,1,*g);
		add_edge(8,12,1,*g);
		add_edge(9,13,1,*g);
		add_edge(9,14,1,*g);
		add_edge(10,11,1,*g);
		add_edge(11,12,1,*g);
		add_edge(13,14,1,*g);
		add_edge(15,16,1,*g);
		add_edge(16,17,1,*g);
		add_edge(17,19,1,*g);
		add_edge(19,20,1,*g);
		add_edge(22,23,1,*g);
		add_edge(18,17,1,*g);
		add_edge(18,19,1,*g);
		add_edge(18,20,1,*g);
		add_edge(18,21,1,*g);
		add_edge(21,22,1,*g);
		add_edge(21,23,1,*g);

		put(poids_sommets,0,3);

		for(int i=1;i<4;i++)
			put(poids_sommets,i,2.5);
		for(int i=4;i<7;i++)
			put(poids_sommets,i,2);
		put(poids_sommets,18,2);
		put(poids_sommets,21,2);
		for(int i=7;i<10;i++)
			put(poids_sommets,i,1.5);
		for(int i=15;i<18;i++)
			put(poids_sommets,i,1.5);
		put(poids_sommets,19,1.5);
		put(poids_sommets,20,1.5);
		put(poids_sommets,22,1.5);
		put(poids_sommets,23,1.5);
		for(int i=10;i<15;i++)
			put(poids_sommets,i,1);

		tie(vertexIt, vertexEnd) = vertices(*g);
		for (; vertexIt != vertexEnd; ++vertexIt)
		{
			cout << *vertexIt << " est connecté avec ";
			tie(neighbourIt, neighbourEnd) = adjacent_vertices(*vertexIt, *g);
			for (; neighbourIt != neighbourEnd; ++neighbourIt)
				cout << *neighbourIt << " ";
			cout<<" et son poids est de "<< poids_sommets[*vertexIt];
			cout << "\n";
		}*/


}

int rand_fini(int a, int b){
	return rand()%(b-a)+a;
}

/**
 * Fonction de recherche d'une valeur dans un tableau.
 * @param tab
 * @param val
 * @return
 */
int recherche_val2(const vector<float> &tab,float val){
	int cpt=0;
	while(tab[cpt]!=val)
		cpt++;
	return cpt;
}

int recherche_val_double(const vector<double> &tab,double val){
	int cpt=0;
	while(tab[cpt]!=val)
		cpt++;
	return cpt;
}

int recherche_val(const Entiers &tab,int val){
	int cpt=0;
	while(tab[cpt]!=val)
		cpt++;
	return cpt;
}
/**
 * @param tab
 * @param i
 * @return
 */
int dichotomie(const Entiers &tab, int i){

	/* déclaration des variables locales à la fonction */
	int id;  //indice de début
	int ifin;  //indice de fin
	int im;  //indice de "milieu"

	/* initialisation de ces variables avant la boucle de recherche */
	id = 0;  //intervalle de recherche compris entre 0...
	ifin = tab.size();  //...et nbVal

	/* boucle de recherche */
	while ((ifin - id) > 1){
		im = (id + ifin)/2;  //on détermine l'indice de milieu
		if(tab[im] > i) ifin = im;  //si la valeur qui est à la case "im" est supérieure à la valeur recherchée, l'indice de fin "ifin" << devient >> l'indice de milieu, ainsi l'intervalle de recherche est restreint lors du prochain tour de boucle
		else id = im;  //sinon l'indice de début << devient >> l'indice de milieu et l'intervalle est de la même façon restreint
	}

	/* test conditionnant la valeur que la fonction va renvoyer */
	if (tab[id] == i) return id;  //si on a trouvé la bonne valeur, on retourne l'indice
	else return -1;  //sinon on retourne -1

}

/**
 * Fonction permettant de supprimer une case d'un tableau.
 * @param tab une référence sur un tableau d'entiers
 * @param i un indice dans tab
 */
void suprim_val(Entiers &tab,int i) {
	tab.erase(tab.begin() + dichotomie(tab,i));
}

/**
 * Détermine si une valeur se trouve dans un tableau.
 * @param tab une référence sur un tableau d'entiers
 * @param val une valeur
 * @return true si la valeur a été trouvée, false sinon
 */
bool In_tab(const Entiers &tab, int val)
{
	for (uint i=0; i < tab.size(); i++) 
		if(tab[i]==val)
			return true;
	return false;
}

bool In_tab_dichotomie(const Entiers &tab, int val)
{
	if(dichotomie(tab,val)!=-1)
		return true;
	else
		return false;
}


void Liste_Voisin(Entiers &P,Entiers &tab,const GraphNonOriente &g)
{

	tie(neighbourIt, neighbourEnd) = adjacent_vertices(P[P.size()-1], g);
	for (; neighbourIt != neighbourEnd; ++neighbourIt)
	{
		if((In_tab(tab,*neighbourIt) == false ) && (In_tab(P,*neighbourIt) == false ))
			tab.push_back(*neighbourIt);
	}
}

int Cout_coupe(Entiers P,int val, GraphNonOriente &g)
{
	int cpt=0;
	P.push_back(val);
	for(int i=0;i<P.size();i++){
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(P[i], g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt){
			if(In_tab(P,*neighbourIt)!=1){
				cpt++;
			}
		}
	}
	return cpt;
}

double Cut_one_cluster(const Entiers &cluster, GraphNonOriente &g, string name)
{
	if(name=="norm"){
		property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),g);
		edge_t e1;
		bool found;
		double cpt=0.;
		for(int i=0;i<cluster.size();i++){
			tie(neighbourIt, neighbourEnd) = adjacent_vertices(cluster[i], g);
			for (; neighbourIt != neighbourEnd; ++neighbourIt){
				tie(e1,found)=edge(vertex(cluster[i],g),vertex(*neighbourIt,g),g);
				if(In_tab(cluster,*neighbourIt)!=1){
					cpt+=get(poids_arc,e1);
				}
			}
		}
		double vol = Cluster_Degree(g,cluster);
		return (cpt/2.)/vol;
	}
	else{
		property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),g);
		edge_t e1;
		bool found;
		double cpt=0.;
		for(int i=0;i<cluster.size();i++){
			tie(neighbourIt, neighbourEnd) = adjacent_vertices(cluster.at(i), g);
			for (; neighbourIt != neighbourEnd; ++neighbourIt){
				tie(e1,found)=edge(vertex(cluster.at(i),g),vertex(*neighbourIt,g),g);
				if(In_tab(cluster,*neighbourIt)!=1){
					cpt+=get(poids_arc,e1);
				}
			}
		}
		return cpt/2.;
	}
}

double Cut_cluster(const EntiersEntiers &tab_cluster,GraphNonOriente &g,string name)
{
	double cpt=0.;
	for(int i=0;i<tab_cluster.size();i++){
		cpt+=Cut_one_cluster(*tab_cluster[i],g,name);
	}
	return cpt;
}

float Cout_coupe_pond(Entiers P,int val, GraphNonOriente &g)
{
	property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),g);
	edge_t e1;
	bool found;
	int cpt=0;

	P.push_back(val);
	for(int i=0;i<P.size();i++){
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(P[i], g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt){
			tie(e1,found)=edge(vertex(P[i],g),vertex(*neighbourIt,g),g);
			if(In_tab(P,*neighbourIt)!=1){
				cpt+=get(poids_arc,e1);
			}
		}
	}
	return cpt;
}


int In_community_dichotomie(const EntiersEntiers &part, int val){
	for (uint i=0; i <part.size() ; i++)
		if (In_tab_dichotomie(*part[i],val))
			return i;
}


/**
 * 
 * @param g
 * @param cluster
 * @return
 */

double Degree(GraphNonOriente &g ,int node){
		property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),g);
		edge_t e1;
		bool found;
		double val=0.;

		tie(neighbourIt, neighbourEnd) = adjacent_vertices(node,g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt){
			tie(e1,found)=edge(vertex(node,g),vertex(*neighbourIt,g),g);
			val+=get(poids_arc,e1);
		}
		return val;
}

double Cluster_Degree(GraphNonOriente &g , const Entiers &cluster){
		property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),g);
		edge_t e1;
		bool found;
		double val=0.;

		for(uint i=0;i<cluster.size();i++){
			val+=Degree(g,cluster.at(i));
		}
		return val;
}

/*double In_modularity(GraphNonOriente &g , const Entiers &cluster){
	property_map<GraphNonOriente,edge_weight_t>::type poids_arc=get(edge_weight_t(),g);
	edge_t e1;
	bool found;
	int val=0;

	for(uint i=0;i<cluster.size();i++){
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(cluster[i],g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt){
			tie(e1,found)=edge(vertex(cluster[i],g),vertex(*neighbourIt,g),g);
			if(In_tab(cluster,*neighbourIt)==1)
				val+=get(poids_arc,e1);
		}
	}
	return val/2.;
}*/

/**
 * 
 * @param g
 * @param cluster
 * @return
 */



/**
 *
 * @param g
 * @param part
 * @return
 */

/*double Modularity(GraphNonOriente &g,const EntiersEntiers &part) {
  double q  = 0.;
  int tmp=num_edges(g);
  for(uint i=0;i<part.size();i++){
	  q+=In_modularity(g,*part[i])/tmp-(Cluster_Degree(g,*part[i])/(2*tmp))*(Cluster_Degree(g,*part[i])/(2*tmp));
  	}

  return q;
}*/

/**
 * 
 * @param part
 * @param val
 * @return
 */



/**
 * Fonction de calcul du gain de modularité de déplacement d'un sommet d'une communoté à une autre !!!!!
 *  ajoute le sommet à part[val] et on calcul la nouvelle modularité
 *  on prend la différence entre la modularité et la nouvouvelle !
 * @param cur_mod
 * @param val
 * @param neight
 * @param node_comm
 * @param part
 * @param g
 */
/*double Modularity_gain(double cur_mod , int val , int neight , int node_comm , EntiersEntiers part , GraphNonOriente &g) {
	double q;
	part[neight]->push_back(val);
	sort(*part[neight]);
	q=Modularity(g,part);

	return q-cur_mod;
}*/

/**
 * Fonction de calcul du gain de modularité de déplacement d'un sommet d'une communoté à une autre !!!!!
 *  ajoute le sommet à part[val] et on calcul la nouvelle modularité
 *  on prend la différence entre la modularité et la nouvouvelle !
 * @param cur_mod
 * @param tmp_community
 * @param neight
 * @param node_comm
 * @param part
 * @param g
 */
/*double Modularity_gain_phase_2(double cur_mod, Entiers tmp_community, int neight, int node_comm, EntiersEntiers part, GraphNonOriente &g) {
	double q;
	for (uint i=0;i<tmp_community.size();i++)
		part[neight]->push_back(tmp_community[i]);
	sort(*part[neight]);
	q = Modularity(g,part);
	return q - cur_mod;
}*/

/**
 * Donne la liste des communautés voisines à celle contenant le sommet val.
 * @param part
 * @param val
 * @param g
 * @return
 */
/*Entiers Neight_community(const EntiersEntiers &part, int val , GraphNonOriente &g){
	Entiers Neight;
	tie(neighbourIt, neighbourEnd) = adjacent_vertices(val, g);
	for (; neighbourIt != neighbourEnd; ++neighbourIt){
		int tmp=In_community(part,*neighbourIt);
		if(In_tab(Neight,tmp)!=1 && In_tab(*part[In_community(part,val)],*neighbourIt)!=1)
			Neight.push_back(tmp);
	}
	sort(Neight);
	return Neight;
}*/

/**
 * 
 * @param part
 * @param community
 * @param g
 * @return
 */
/*Entiers Part_Neight_community(const EntiersEntiers &part,int community, GraphNonOriente &g){
	Entiers Neight;
	for(uint i =0;i<part[community]->size();i++){
		tie(neighbourIt, neighbourEnd) = adjacent_vertices(part[community]->at(i), g);
		for (; neighbourIt != neighbourEnd; ++neighbourIt){
			int tmp=In_community(part,*neighbourIt);
			if(In_tab(Neight,tmp)!=1 && tmp!=community)
				Neight.push_back(tmp);
		}
	}
	sort(Neight);
	return Neight;
}*/
