#include "lib_utils.h"


/**
 * \brief Scales the data X into the unit hypercube
 *
 * sometimes data has a funny scale and it is outside of the libigl-viewer-frustrum
 * this is a convenience function to scale the data into the unit hypercube
 *
 * \params[in/out] X The data points organized as a n x d matrix
 */
void scale_to_unit_cube(Eigen::MatrixXd &X, double &s)
{
    if (s == 0)
    {
        s = 1.0 / X.colwise().maxCoeff().maxCoeff();
    }

    X = (X.array() * s).matrix();
}

/**
 * \brief Finds the length of the smallest edge in the mesh
 *
 * why? the convergence properties of the Laplacian depend on the smallest edge length
 *  so in case you need that you can scale your mesh accordingly
 *
 * \params[in] points The data points organized as a n x d matrix
 * \params[in] faces The faces (connectivity) organized as a n x d matrix
 */
double find_min_edge_length(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces)
{
    double min = std::numeric_limits<double>::max();

    for(int face_idx=0; face_idx < faces.rows(); face_idx++)
    {
        int i = faces(face_idx, 0);
        int j = faces(face_idx, 1);
        int k = faces(face_idx, 2);

        Eigen::Vector3d v_i = points.row(i);
        Eigen::Vector3d v_j = points.row(j);
        Eigen::Vector3d v_k = points.row(k);

        Eigen::Vector3d e_ij = v_j - v_i;
        Eigen::Vector3d e_jk = v_k - v_j;
        Eigen::Vector3d e_ki = v_i - v_k;

        min = std::min(e_ij.norm(), min);
        min = std::min(e_jk.norm(), min);
        min = std::min(e_ki.norm(), min);
    }

    return min;
}