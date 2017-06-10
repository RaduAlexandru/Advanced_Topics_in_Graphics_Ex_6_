#ifndef ATCG1_LIB_UTILS_H
#define ATCG1_LIB_UTILS_H

#include <Eigen/Core>

void scale_to_unit_cube(Eigen::MatrixXd &X, double &s);

double find_min_edge_length(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces);

#endif //ATCG1_LIB_UTILS_H

