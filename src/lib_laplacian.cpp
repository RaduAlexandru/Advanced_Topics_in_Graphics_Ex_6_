#include <iostream>
#include "lib_laplacian.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Geometry>

/**
 * \brief Build the graph laplacian (also called umbrella operator) as a sparse matrix
 *
 * Ref:
 *  -Desbrun et al.: "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow"
 *  -(Equation 7)
 *
 * what it does:
 *  -construct the graph laplacian from the connectivity encoded in the faces
 *
 * TODO:
 *  -construct the sparse symmetric adjacency matrix A
 *  -compute the diagonal matrix D with weights d_ii = sum_j a_ij
 *  -construct the graph laplacian from A and D
 *
 *  \param[in] faces The triangles encoding the connectivity
 *  \param[in/out] M The lumped mass matrix
 *  \param[in/out] S The stiffness matrix
 */
Eigen::SparseMatrix<double> build_graph_laplacian(const Eigen::MatrixXi &faces, Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &S)
{
    int nv = faces.maxCoeff() + 1;
    Eigen::SparseMatrix<double> L(nv, nv);
    M.resize(nv, nv);
    S.resize(nv, nv);




    return L;
}

/**
 * \brief Safe version for the acos function
 *
 *  -You should use this function!
 *
 *  \param[in] value The parameter in -1.0 to 1.0
 *  \param[out] acos(value) The computed acos value
 */
inline double safe_acos(const double value)
{
    if(value > 1.0)
    {
        std::cout <<" >1" << std::endl;
        return acos(1.0);
    }
    else if(value < -1.0)
    {
        std::cout <<" <-1.0" << std::endl;
        return acos(-1.0);
    }
    else
    {
        return acos(value);
    }
}


/**
 * \brief Computes a numerically stable version of herons formula for the area of a triangle from the metric
 *
 *  -You should use this function!
 *
 *  \param[in] a,b,c The edge lengths in unsorted order
 *  \param[out] area The stable area for spikey triangles
 */
inline double triangle_area_from_metric(double a, double b, double c)
{
    /*
     * numerically stable version of herons formula
    */

    double max = std::max(a, std::max(b, c));
    double min = std::min(a, std::min(b, c));
    double last = a + b + c - max - min;

    a = max; b = last; c = min;
    assert((a >= b) && (b >= c)); // ordering of the edges by length

    double term1 = a + (b + c);
    double term2 = c - (a - b);
    double term3 = c + (a - b);
    double term4 = a + (b - c);

    double area = std::sqrt(term1 * term2 * term3 * term4) / 4.0;

    if(area < 1e-8)
    {
        area = 1e-8;
        std::cout << "small area < 1e-8!" << std::endl;
    }
    return area;
}

/**
 * \brief Computes a numerically stable version of law of cosines for the angle from the metric
 *
 *  -You should use this function!
 *
 *  \param[in] a,b,c The edge lengths in unsorted order
 *  \param[out] alpha The angle between a and b, opposite to c
 */
inline double angle_from_metric(double a, double b, double c)
{
    /* numerically stable version of law of cosines
     * angle between a and b, opposite to c
    */

    double alpha = safe_acos((a*a + b*b - c*c) / (2.0 * a * b));

    if (alpha < 1e-8)
    {
        alpha = std::sqrt((c*c - (a - b)*(a - b)) / (2.0 * a * b));
        std::cout << "small angle < 1e-8!" << std::endl;
    }
    return alpha;
}

/**
 * \brief Build the cotan laplacian as a sparse matrix
 *
 * Ref:
 *  -Desbrun et al.: "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow"
 *  -(Equation 14)
 *
 * what it does:
 *  -construct the cotan laplacian from the connectivity encoded in the faces
 *  -and the specific cotan-weights for each edge
 *  -the normalize parameter controls wether you want the normalized laplacian or weight each row by the inverse of the vertex area measure
 *
 * TODO:
 *  -iterate over each face
 *  -compute the angle at each vertex
 *  -compute the triangle area
 *  -compute the cotan weight from the angles
 *  -fill a sparse matrix L with the cotan weights by using Eigen::Triplet<double>
 *  -fill a diagonal matrix D with the inverse area weights
 *
 *  -make sure that W is symmetric!
 *   (pseudo code:)
 *  -then compute L = W - M
 *  -then normalize L = M^(-1) * L with
 *      M = rowsum or with M = areaweights
 *
 *  \param[in] faces The triangles encoding the connectivity
 *  \param[in/out] M The lumped mass matrix
 *  \param[in/out] S The stiffness matrix
 */
Eigen::SparseMatrix<double> build_cotan_laplacian(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, Eigen::SparseMatrix<double> &M, Eigen::SparseMatrix<double> &S, const bool normalize)
{
    int nv = points.rows();
    Eigen::SparseMatrix<double> L(nv, nv);
    M.resize(nv, nv);
    S.resize(nv, nv);
    return L;
}

/**
 * \brief Explicit euler steps for smoothing
 *
 * Ref:
 *  -Desbrun et al.: "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow"
 *
 * TODO:
 *  -iterate the smoothing nsteps times
 *  -dont forget to compute the steplength!
 *  -if you compute the smoothing with the cotan laplacian, remember to recompute the cotan laplacian after
 *      the iteration because the metric changed!
 *      Otherwise your flow is biased and not a true mean curvature flow!
 *
 *  \param[in] points The points of the mesh (n x 3)
 *  \param[in] faces The faces of the mesh (m x 3)
 *  \param[in] nsteps The number of steps to iterate
 *  \param[in] time The time value for the diffusion process
 *  \param[in] L The Laplacian
 *  \param[in] recompute Wether we have to recompute the cotan Laplacian
 *  \param[in] normalize Wether we have to normalize the cotan Laplacian
 *  \param[out] smooth_points The result of the diffusion process
 */
Eigen::MatrixXd smooth_explicit(const Eigen::MatrixXd points,
                                const Eigen::MatrixXi &faces,
                                const int nsteps,
                                const double time,
                                Eigen::SparseMatrix<double> &L,
                                const bool recompute,
                                const bool normalize)
{
    Eigen::MatrixXd smooth_points(points.rows(), points.cols());
    return smooth_points;
}

/**
 * \brief Implicit euler steps for smoothing
 *
 * Ref:
 *  -Desbrun et al.: "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow"
 *
 * TODO:
 *  -iterate the smoothing nsteps times
 *  -dont forget to compute the steplength!
 *  -if you compute the smoothing with the cotan laplacian, remember to recompute the cotan laplacian after
 *      the iteration because the metric changed!
 *      Otherwise your flow is biased and not a true mean curvature flow!
 *
 *  -You can use Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver; to solve the equation!
 *16
 *  \param[in] points The points of the mesh (n x 3)
 *  \param[in] faces The faces of the mesh (m x 3)
 *  \param[in] nsteps The number of steps to iterate
 *  \param[in] time The time value for the diffusion process
 *  \param[in] L The Laplacian
 *  \param[in] recompute Wether we have to recompute the cotan Laplacian
 *  \param[in] normalize Wether we have to normalize the cotan Laplacian
 *  \param[out] smooth_points The result of the diffusion process
 */
Eigen::MatrixXd smooth_implicit(const Eigen::MatrixXd points, const Eigen::MatrixXi &faces, const int nsteps, const double time, Eigen::SparseMatrix<double> &L, const bool recompute, const bool normalize)
{
    Eigen::MatrixXd smooth_points(points.rows(), points.cols());
    return smooth_points;
}
