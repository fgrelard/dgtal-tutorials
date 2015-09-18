#ifndef __DIANA__
#define __DIANA__
#include <vector>
using namespace std;

namespace Agnes {
	template <typename Point>
	double euclideanDistance(const Point& first, const Point& second);

	template <typename Point>
	int* diameter(const vector<vector<Point> >& clusters);

	template <typename Point>
	void clusterBelonging(int*, vector<vector<Point> >& clusters);

	template <typename Point>
	void mainLoop(vector<vector<Point > >& clusters);
}

#include "clustering/diana.hpp"
#include <math.h>
#include <iostream>
#include <algorithm>
#include <limits>
using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////

template <typename Point>
double Agnes::euclideanDistance(const Point& x, const Point& y) {
	double dist = sqrt(pow(x[0]-y[0], 2) + pow(x[1]-y[1], 2) + pow(x[2]-y[2], 2));
	return dist;
}

template <typename Point>
int* Agnes::diameter(const vector<vector<Point> >& cluster) {
	double min = numeric_limits<double>::max();
	double diameter = 0.0;
	int* indice_diameter = NULL;
	vector<Point> mean_vector;
	for (int i = 0; i < cluster.size(); i++) {
		Point point(min, min, min);
		for (int j = 0; j < cluster[i].size(); j++) {
			if (cluster[i][j] < point)
				point = cluster[i][j];
		}
		mean_vector.push_back(point);
	} 

	if (mean_vector.size() == 0)
		return indice_diameter;
	//critere d'arret : % distance maximale
	double seuil = euclideanDistance(*max_element(mean_vector.begin(), mean_vector.end()), *min_element(mean_vector.begin(), mean_vector.end()))+1;
	for (int i = 0; i<mean_vector.size(); i++){
		for (int j = i+1; j<mean_vector.size(); j++){
			diameter = euclideanDistance(mean_vector[i],mean_vector[j]);
			if (diameter<min && diameter < seuil){
				if (!indice_diameter)
					indice_diameter = new int[2];
				min = diameter;
				indice_diameter[0] = i;
				indice_diameter[1] = j;
			}
		}
	}
	//les indices des clusters a fusionner sont retournés
	return indice_diameter;
}


template <typename Point>
void Agnes::clusterBelonging(int* indice, vector<vector<Point> >& clusters) {
	// correspond au nombre de proteines
	vector<Point> nouveau_cluster;

	for (int i = 0; i <clusters[indice[0]].size(); i++) {
		nouveau_cluster.push_back(clusters[indice[0]][i]);
	}
	for (int i = 0; i <clusters[indice[1]].size(); i++) {
		nouveau_cluster.push_back(clusters[indice[1]][i]);
	}
	clusters.erase(find(clusters.begin(), clusters.end(), clusters[indice[0]]));
	clusters.erase(find(clusters.begin(), clusters.end(), clusters[indice[1]-1]));
	clusters.push_back(nouveau_cluster);


	double sum_splinter = 0.0;
	double sum_remain = 0.0;
	double size = clusters.size();
  
	//Tous les clusters sont parcourus (sauf le nouveau)
	for (int i = 0; i < size - 1; i++) {
		
		for (int j = 0; j < clusters[i].size(); j++) {
			sum_splinter = 0.0;
			sum_remain = 0.0;
			
			for (int k = 0; k < clusters[i].size(); k++) {
				if (j != k) {
					sum_remain += euclideanDistance(clusters[i][j], clusters[i][k]);
				}
			}
			for (int k = 0; k < clusters[size - 1].size(); k++) {
				sum_splinter += euclideanDistance(clusters[i][j], clusters[size - 1][k]);
			}
			sum_splinter /= clusters[size - 1].size();
			sum_remain /= clusters[i].size();

			//Si la moyenne des distances du point considéré est plus faible dans le nouveau, on l'ajoute
			if (sum_splinter < sum_remain) {
				clusters[size - 1].push_back(clusters[i][j]);
				clusters[i].erase(clusters[i].begin()+j);
			}
		}
	}
}

template <typename Point>
void Agnes::mainLoop(vector<vector<Point> >& clusters) {
	int* indice_diameter=diameter(clusters);
	while (indice_diameter){  
		clusterBelonging(indice_diameter, clusters);
		indice_diameter=diameter(clusters);
	}
}


#endif
