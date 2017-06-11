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


    // //Test
    // nv=4;
    // Eigen::MatrixXi faces_test(2,3);
    // faces_test(0,0)=0;
    // faces_test(0,1)=1;
    // faces_test(0,2)=2;
    // faces_test(1,0)=1;
    // faces_test(1,1)=3;
    // faces_test(1,2)=2;
    // std::cout << "face test is \n" << faces_test << '\n';

    //adjacency matrix
    Eigen::SparseMatrix<double> adjacency_matrix(nv, nv);
    std::vector<Eigen::Triplet<double>> adjacency_coeff;
    // Eigen::MatrixXd dense_adjacency= Eigen::MatrixXd::Zero(nv,nv);
    for (size_t i = 0; i < faces.rows(); i++) {
      //there are 3 vertices to a triangle and we need the adyancency between all of them
      // vertces 0 and 1
      // dense_adjacency(faces(i,0), faces(i,1)) = 1;
      // dense_adjacency(faces(i,1), faces(i,0)) = 1;
      //
      // // //vertex 1 and 2
      // dense_adjacency(faces(i,1), faces(i,2)) = 1;
      // dense_adjacency(faces(i,2), faces(i,1)) = 1;
      //
      // // //vertex 2 and 0
      // dense_adjacency(faces(i,2), faces(i,0)) = 1;
      // dense_adjacency(faces(i,0), faces(i,2)) = 1;

      // // vertces 0 and 1
      // if (adjacency_matrix.coeff(faces_test(i,0), faces_test(i,1))!=1 ){
      //    adjacency_matrix.coeffRef(faces_test(i,0), faces_test(i,1))= 1;
      // }
      // if (adjacency_matrix.coeff(faces_test(i,1), faces_test(i,0))!=1 ){
      //    adjacency_matrix.coeffRef(faces_test(i,1), faces_test(i,0))= 1;
      // }
      //
      // // // //vertex 1 and 2
      // if (adjacency_matrix.coeff(faces_test(i,1), faces_test(i,2))!=1 ){
      //    adjacency_matrix.coeffRef(faces_test(i,1), faces_test(i,2))= 1;
      // }
      // if (adjacency_matrix.coeff(faces_test(i,2), faces_test(i,1))!=1 ){
      //    adjacency_matrix.coeffRef(faces_test(i,2), faces_test(i,1))= 1;
      // }
      //
      // //vertex 2 and 0
      // if (adjacency_matrix.coeff(faces_test(i,2), faces_test(i,0))!=1 ){
      //    adjacency_matrix.coeffRef(faces_test(i,2), faces_test(i,0))= 1;
      // }
      // if (adjacency_matrix.coeff(faces_test(i,0), faces_test(i,2))!=1 ){
      //    adjacency_matrix.coeffRef(faces_test(i,0), faces_test(i,2))= 1;
      // }

      // vertex 0 and 1
      adjacency_coeff.emplace_back(faces(i,0), faces(i,1), 1);
      adjacency_coeff.emplace_back(faces(i,1), faces(i,0), 1);

      // vertex 1 and 2
      adjacency_coeff.emplace_back(faces(i,1), faces(i,2), 1);
      adjacency_coeff.emplace_back(faces(i,2), faces(i,1), 1);

      // vertex 2 and 0
      adjacency_coeff.emplace_back(faces(i,2), faces(i,0), 1);
      adjacency_coeff.emplace_back(faces(i,0), faces(i,2), 1);

    }
    adjacency_matrix.setFromTriplets(adjacency_coeff.begin(), adjacency_coeff.end());

    //Set the adjacency matrix to 1
    for (size_t i = 0; i < adjacency_matrix.rows(); i++) {
      for (size_t j = 0; j < adjacency_matrix.cols(); j++) {
        if (adjacency_matrix.coeff(i,j)>1) {
          adjacency_matrix.coeffRef(i,j)=1;
        }
      }
    }




    // Eigen::MatrixXd dense_adjacency = Eigen::MatrixXd(adjacency_matrix);
    // std::cout << "adyancecy matrix is \n " << dense_adjacency << '\n';
    // adjacency_matrix=dense_adjacency.sparseView();

    // // //Degree matrix
    // Eigen::SparseMatrix<double> degree_matrix(nv, nv);
    // std::vector<Eigen::Triplet<double>> degree_coeff;
    // for (size_t i = 0; i < dense_adjacency.rows(); i++) {
    //   int degree=0;
    //   for (size_t j = 0; j < dense_adjacency.cols(); j++) {
    //     if (dense_adjacency(i,j)!=0 ){
    //       degree++;
    //     }
    //   }
    //   degree_coeff.emplace_back(i, i, degree);
    // }
    // degree_matrix.setFromTriplets(degree_coeff.begin(), degree_coeff.end());
    // // Eigen::MatrixXd dense_degree = Eigen::MatrixXd(degree_matrix);
    // // std::cout << "degree_matrix is \n " << dense_degree << '\n';
    //
    // L=0.5*(adjacency_matrix-degree_matrix);


    //degree matrix
    Eigen::SparseMatrix<double> degree_matrix(nv, nv);
    std::vector<Eigen::Triplet<double>> degree_coeff;
    for (size_t i = 0; i < adjacency_matrix.rows(); i++) {
      int degree=0;
      for (size_t j = 0; j < adjacency_matrix.cols(); j++) {
        if (adjacency_matrix.coeff(i,j)!=0 ){
          degree++;
        }
      }
      degree_coeff.emplace_back(i, i, degree);
    }
    degree_matrix.setFromTriplets(degree_coeff.begin(), degree_coeff.end());
    // Eigen::MatrixXd dense_degree = Eigen::MatrixXd(degree_matrix);
    // std::cout << "degree_matrix is \n " << dense_degree << '\n';

    L=0.5*(adjacency_matrix-degree_matrix);

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

    // Create the adjacency matrix
    Eigen::SparseMatrix<double> adjacency_matrix(nv, nv);
    for (size_t i = 0; i < faces.rows(); i++)
    {
      int x = faces(i,0);
      int y = faces(i,1);
      int z = faces(i,2);
      adjacency_matrix.coeffRef(x,y) = 1;
      adjacency_matrix.coeffRef(y,x) = 1;
      adjacency_matrix.coeffRef(x,z) = 1;
      adjacency_matrix.coeffRef(z,x) = 1;
      adjacency_matrix.coeffRef(y,z) = 1;
      adjacency_matrix.coeffRef(z,y) = 1;
}



    for (size_t i = 0; i < faces.rows(); i++) {
      Eigen::VectorXd p1=points.row(faces(i,0));
      Eigen::VectorXd p2=points.row(faces(i,1));
      Eigen::VectorXd p3=points.row(faces(i,2));

      std::vector<float> angles(3);

      //angle at p1
      Eigen::VectorXd p1_l1= (p2-p1).normalized();
      Eigen::VectorXd p1_l2= (p3-p1).normalized();
      angles[0]= safe_acos(p1_l1.dot(p1_l2));

      Eigen::VectorXd p2_l1= (p1-p2).normalized();
      Eigen::VectorXd p2_l2= (p3-p2).normalized();
      angles[1]= safe_acos(p2_l1.dot(p2_l2));

      Eigen::VectorXd p3_l1= (p2-p3).normalized();
      Eigen::VectorXd p3_l2= (p1-p3).normalized();
      angles[2]= safe_acos(p3_l1.dot(p3_l2));

      // float sum_of_elems=0.0;
      // for (auto& n : angles)
      //   sum_of_elems += n;
      // std::cout << "angles sum is" <<sum_of_elems << '\n';

      //Triangle area
      float l1_norm= (p1-p2).norm();
      float l2_norm= (p1-p3).norm();
      float l3_norm= (p2-p3).norm();
      float area= triangle_area_from_metric(l1_norm,l2_norm,l3_norm);

      // std::cout << "angle is" << angles[0] << '\n';

    }



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

    smooth_points=points;
    for (size_t i = 0; i < 50; i++) {
      smooth_points= smooth_points+ 0.25*L*smooth_points;
    }



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
    smooth_points=points;

    float steplength=0.25;

    //Ax=b
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;

    //A
    Eigen::SparseMatrix<double> I (points.rows(), points.rows());
    I.setIdentity();
    Eigen::SparseMatrix<double> A = I-steplength*L;



    //solve
    solver.compute(A);
    for (size_t i = 0; i < 50; i++) {
      smooth_points = solver.solve(smooth_points);
    }





    return smooth_points;
}
