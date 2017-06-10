#include "lib_knn_query.h"
#include <nanoflann.hpp>
#include <cstdint>

/**
 * \brief Find the k nearest neighbors for each point in query_points in target_points using a KDTree
 *
 * what it does:
 *  -build a kd_tree for the target_points
 *  -find for each point in query_points the k closest points in target_points
 *  -put the found indices and distances into the out_indices and out_dists matrices
 *
 * \param[in] target_points The target points organized as a n x d matrix
 * \param[in] query_points The query points organized as a m x d matrix
 * \param[in] k The number of nearest neighbors
 * \param[out] out_indices The found indices in the target points
 * \param[out] out_dists The distances to the found points in the target points
 */

void findKNN(const Eigen::MatrixXd &target_points,
             const Eigen::MatrixXd &query_points,
             const int k,
             Eigen::MatrixXi &out_indices,
             Eigen::MatrixXd &out_dists)
{
    typedef nanoflann::KDTreeEigenMatrixAdaptor< Eigen::MatrixXd >  my_kd_tree_t;

    my_kd_tree_t mat_index(target_points.cols(), target_points, 10 /*max_leaf*/);
    mat_index.index->buildIndex();

    std::vector<int64_t> tmp_indices(k);
    std::vector<double> tmp_dists(k);
    std::vector<double> query_pt(3);

    for(int i=0; i < query_points.rows(); i++){
        query_pt[0] = query_points.row(i)[0];
        query_pt[1] = query_points.row(i)[1];
        query_pt[2] = query_points.row(i)[2];

        mat_index.index->knnSearch(&query_pt[0], k, &tmp_indices[0], &tmp_dists[0]);

        for(int j=0;j < k; j++)
        {
            out_dists(i,j) = tmp_dists[j];
            out_indices(i,j) = tmp_indices[j];
        }
    }
}