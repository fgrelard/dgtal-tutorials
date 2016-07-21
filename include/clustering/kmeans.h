#ifndef KMEANS_H
#define KMEANS_H

#include <vector>
#include <limits>
#include "Statistics.h"

namespace KMeans {

        template <typename Scalar>
        std::vector<std::vector<Scalar > > kmeansAlgorithm(const std::vector<Scalar>& initialData, int k);

        template <typename Scalar>
        std::vector<std::vector<Scalar> > initialize2Clusters(const std::vector<Scalar>& initialData);

        template <typename Scalar>
        void assignClasses(std::vector<std::vector<Scalar> > & clusters, const std::vector<Scalar>& initialData, const std::vector<Scalar>& centroids);

        template <typename Scalar>
        void recomputeCentroids( std::vector<Scalar>& centroids, const std::vector<std::vector<Scalar> > & clusters);
};

template <typename Scalar>
std::vector<std::vector<Scalar > > KMeans::kmeansAlgorithm(const std::vector<Scalar>& initialData, int  k) {
        std::vector<std::vector<Scalar> > clusters = initialize2Clusters(initialData);

        std::vector<Scalar> centroids;
        for (const std::vector<Scalar>& point : clusters) {
                centroids.push_back(point[0]);
        }

        int cpt = 0;
        bool convergence = false;
        while (!convergence && cpt < 1000) {
                std::vector<Scalar> previousCentroids = centroids;

                assignClasses(clusters, initialData, centroids);
                recomputeCentroids(centroids, clusters);
                cpt++;

                if (centroids == previousCentroids)
                        convergence = true;
        }
        return clusters;
}

template <typename Scalar>
std::vector<std::vector<Scalar> > KMeans::initialize2Clusters(const std::vector<Scalar>& initialData) {
        std::vector<std::vector<Scalar> > clusters;
        Scalar minVal = *min_element(initialData.begin(), initialData.end());
        Scalar maxVal = *max_element(initialData.begin(), initialData.end());
        clusters.push_back({minVal});
        clusters.push_back({maxVal});
        return clusters;
}

template <typename Scalar>
void KMeans::assignClasses(std::vector<std::vector<Scalar> > & clusters, const std::vector<Scalar>& initialData, const std::vector<Scalar>& centroids) {
        int  k = clusters.size();
        for (int i  = 0; i < k; i++) {
                clusters[i].clear();
        }

        for (Scalar data : initialData) {
                Scalar distance = std::numeric_limits<Scalar>::max();
                int index = 0;
                for (int i = 0; i < k; i++) {
                        Scalar currentDistance = sqrt(pow((data - centroids[i]), 2));
                        if (currentDistance < distance) {
                                distance = currentDistance;
                                index = i;
                        }
                }
                clusters[index].push_back(data);
        }
}

 template <typename Scalar>
 void KMeans::recomputeCentroids( std::vector<Scalar>& centroids, const std::vector<std::vector<Scalar> > & clusters) {
         for (int i = 0, end = clusters.size(); i < end; i++) {
                 double mean = Statistics::mean(clusters[i]);
                 centroids[i] = mean;
         }
 }

#endif
