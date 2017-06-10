#ifndef ATCG1_LIB_LAPLACIAN_H
#define ATCG1_LIB_LAPLACIAN_H

#include <Eigen/Sparse>

Eigen::SparseMatrix<double> build_graph_laplacian(const Eigen::MatrixXi &faces, Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &S);

Eigen::SparseMatrix<double> build_cotan_laplacian(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &S, const bool normalize);

Eigen::MatrixXd smooth_explicit(const Eigen::MatrixXd points, const Eigen::MatrixXi &faces, const int nsteps, const double time, Eigen::SparseMatrix<double> &L, const bool recompute=false, const bool normalize=false);

Eigen::MatrixXd smooth_implicit(const Eigen::MatrixXd points, const Eigen::MatrixXi &faces, const int nsteps, const double time, Eigen::SparseMatrix<double> &L, const bool recompute=false, const bool normalize=false);

#endif //ATCG1_LIB_LAPLACIAN_H
