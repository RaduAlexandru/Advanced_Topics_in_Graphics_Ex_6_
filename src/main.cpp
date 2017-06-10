#include <cstdlib>
#include <iostream>

//set this to supress libigl viewer help
#define IGL_VIEWER_VIEWER_QUIET

#include <igl/viewer/Viewer.h>
#include <igl/readOBJ.h>
#include <igl/parula.h>

#include <igl/boundary_facets.h>
#include <igl/colon.h>
#include <igl/unique.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/setdiff.h>

#include <Eigen/Core>
#include "lib_laplacian.h"
#include "lib_utils.h"

/**
 * \brief First exercise
 *
 * what it does:
 *  -load a mesh
 *  -complete the build_graph_laplacian function in src
 *  -complete the build_cotan_laplacian function (with and without normalization) in src
 *  -smooth the surface with each laplacian for 3 explicit steps with time 400.0 (this should usually break down)
 *  -smooth the surface with each laplacian for 3 implicit steps with time 400.0 (this should give you nice smooth surfaces)
 *
 *  -plot the results for the different smoothing parameters
 *
 *  \param[in] filename1 The path to the mesh
 */
void exercise1(const std::string &filename1){
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    bool load_success_target = igl::readOBJ(filename1, V, F);
    if (!load_success_target)
    {
        std::cerr << "could not load file: " << filename1 << std::endl;
        return;
    }

    double min_edge_length = find_min_edge_length(V, F);
    double s = 1.0 / min_edge_length;
    scale_to_unit_cube(V, s);

    Eigen::SparseMatrix<double> M_graph, S_graph;
    Eigen::SparseMatrix<double> L_graph = build_graph_laplacian(F, M_graph, S_graph);

    Eigen::SparseMatrix<double> M_cot, S_cot;
    Eigen::SparseMatrix<double> L_cot = build_cotan_laplacian(V, F, M_cot, S_cot, false);

    Eigen::SparseMatrix<double> M_norm, S_norm;
    Eigen::SparseMatrix<double> L_norm = build_cotan_laplacian(V, F, M_norm, S_norm, true);

    igl::viewer::Viewer viewer;
    std::cout << "graph laplace, explicit smoothing, 3 steps, time: 400.0" << std::endl;
    Eigen::MatrixXd smooth_expl_graph_V = smooth_explicit(V, F, 3, 400.0, L_graph, false);
    viewer.data.set_mesh(smooth_expl_graph_V, F);
    viewer.launch();

    viewer.data.clear();
    std::cout << "graph laplace, implicit smoothing, 3 steps, time: 400.0" << std::endl;
    Eigen::MatrixXd smooth_impl_graph_V = smooth_implicit(V, F, 3, 400.0, L_graph, false);
    viewer.data.set_mesh(smooth_impl_graph_V, F);
    viewer.launch();

    viewer.data.clear();
    std::cout << "cotan laplace, explicit smoothing, 3 steps, time: 400.0" << std::endl;
    Eigen::MatrixXd smooth_expl_cotan_V = smooth_explicit(V, F, 3, 400.0, L_cot, true, false);
    viewer.data.set_mesh(smooth_expl_cotan_V, F);
    viewer.launch();

    viewer.data.clear();
    std::cout << "cotan laplace, implicit smoothing, 3 steps, time: 400.0" << std::endl;
    Eigen::MatrixXd smooth_impl_cotan_V = smooth_implicit(V, F, 3, 400.0, L_cot, true, false);
    viewer.data.set_mesh(smooth_impl_cotan_V, F);
    viewer.launch();

    viewer.data.clear();
    std::cout << "normalized cotan laplace, explicit smoothing, 3 steps, time: 400.0" << std::endl;
    Eigen::MatrixXd smooth_expl_cotan_V_norm = smooth_explicit(V, F, 3, 400.0, L_norm, true, true);
    viewer.data.set_mesh(smooth_expl_cotan_V_norm, F);
    viewer.launch();

    viewer.data.clear();
    std::cout << "normalized cotan laplace, implicit smoothing, 3 steps, time: 400.0" << std::endl;
    Eigen::MatrixXd smooth_impl_cotan_V_norm = smooth_implicit(V, F, 3, 400.0, L_norm, true, true);
    viewer.data.set_mesh(smooth_impl_cotan_V_norm, F);
    viewer.launch();
}

/**
 * \brief Computes the mean curvature by applying the Laplacian to the vertices
 *
 * TODO:
 *  -apply a Laplacian matrix to the points
 *  -compute the mean curvature from the norm of the delta coordinates
 *
 *  \param[in] points The points (n x 3)
 *  \param[in] L The Laplacian
 */
Eigen::VectorXd compute_mean_curvature(const Eigen::MatrixXd &points, const Eigen::SparseMatrix<double> &L)
{
    Eigen::VectorXd mean_curv(points.rows());
    return mean_curv;
}

/**
 * \brief Second exercise
 *
 * what it does:
 *  -load a mesh
 *  -Show the mean curvature
 *
 *  \param[in] filename1 The path to the mesh
 */
void exercise2(const std::string &filename)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    bool load_success_target = igl::readOBJ(filename, V, F);
    if (!load_success_target) {
        std::cerr << "could not load file: " << filename << std::endl;
        return;
    }

    double min_edge_length = find_min_edge_length(V, F);
    double s = 1.0 / min_edge_length;
    scale_to_unit_cube(V, s);



    igl::viewer::Viewer viewer;

    Eigen::SparseMatrix<double> M_graph, S_graph;
    Eigen::SparseMatrix<double> L_graph = build_graph_laplacian(F, M_graph, S_graph);

    Eigen::SparseMatrix<double> M_cot, S_cot;
    Eigen::SparseMatrix<double> L_cot = build_cotan_laplacian(V, F, M_cot, S_cot, false);

    Eigen::SparseMatrix<double> M_norm, S_norm;
    Eigen::SparseMatrix<double> L_norm = build_cotan_laplacian(V, F, M_norm, S_norm, true);

    Eigen::MatrixXd colors;
    Eigen::VectorXd mean_curv_graph = compute_mean_curvature(V, L_graph);
    igl::parula(mean_curv_graph, true, colors);

    std::cout << "graph laplace, mean curvature" << std::endl;
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(colors);
    viewer.launch();

    viewer.data.clear();

    Eigen::VectorXd mean_curv_cot = compute_mean_curvature(V, L_cot);
    igl::parula(mean_curv_cot, true, colors);

    std::cout << "cotan laplace, mean curvature" << std::endl;
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(colors);
    viewer.launch();

    viewer.data.clear();

    Eigen::VectorXd mean_curv_norm = compute_mean_curvature(V, L_norm);
    igl::parula(mean_curv_norm, true, colors);

    std::cout << "cotan normalized laplace, mean curvature" << std::endl;
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(colors);
    viewer.launch();
}

/**
 * \brief detect the boundary vertices for a given mesh
 *
 * TODO:
 *  -compute all the indices at the boundary of a mesh
 *  -make them unique
 *  -compute the boundary indices from the edges
 *  -compute the interior indices with igl::setdiff, you can use igl::colon to generate a list of all the indices
 *
 *  -plot the results for the different smoothing parameters
 *
 *  \param[in] points The points of a mesh (n x 3)
 *  \param[in] faces The faces of a mesh (m x 3)
 *  \param[in/out] boundary_indices The indices of the boundary (as vector)
 *  \param[in/out] interior_indices The indices of the interior (as vector)
 */
void detect_boundary_vertices(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, Eigen::VectorXi &boundary_indices, Eigen::VectorXi &interior_indices)
{
}

/**
 * \brief Increase the boundary by region growing
 *
 * TODO:
 *  -increase the k-ring of boundary vertices to the k+1 ring of the boundary vertices
 *  -region growing can be done by iteratively spreading values from selected pooints to neighboring points with the laplacian
 *
 *  \param[in] points The points of a mesh (n x 3)
 *  \param[in] L The Laplacian
 *  \param[in/out] boundary_indices The new indices of the increased set of the boundary (as vector)
 *  \param[in/out] interior_indices The new indices of the decreased set of the interior (as vector)
 */
void increase_ring_boundary(const Eigen::MatrixXd &points, const Eigen::SparseMatrix<double> &L,
                            Eigen::VectorXi &boundary_indices, Eigen::VectorXi &interior_indices)
{
}

/**
 * \brief Constraint fairing by solving a linear system
 *
 * TODO:
 *  -compute the smoothness degree k at the boundary by computing a higher order laplacian (k=1, 2, 3)
 *  -construct a linear system by reordering the entries of the laplacian as:
 *      _     _
 *     |A    B|
 *     |------|
 *     |0 | Id|
 *     |_    _|
 *
 *     where A is the laplacian part for the interior vertices
 *     and B is the laplacian part for the interior vertices connected to the boundary vertices
 *
 *  -You can compute them by using igl::slice and compose the linear system with igl::slice_into
 *  -Hint: There is some index juggling involved here!
 *
 *  \param[in] L The Laplacian
 *  \param[in/out] boundary_indices The new indices of the increased set of the boundary (as vector)
 *  \param[in/out] interior_indices The new indices of the decreased set of the interior (as vector)
 *  \param[in] k The smoothness degree
 */
Eigen::SparseMatrix<double> construct_laplacian_fairing_system(const Eigen::SparseMatrix<double> &L,
                                                               const Eigen::VectorXi &boundary_indices,
                                                               const Eigen::VectorXi &interior_indices, int k = 1)
{
    Eigen::SparseMatrix<double> system(L.rows(), L.cols());
    return system;
}


/**
 * \brief Fairing function to calculated the fair surface connecting the fixed boundary
 *
 * TODO:
 *  -complete the constraint_fairing function which will give you the operator of your linear system
 *  -construct the righthandside by setting all points in the interior to 0 and fix the boundary at the current positions
 *  -you can solve the laplace equation with     Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver(op);
 *      which was a stable solver in the test runs.
 *
 *  -Note that for k=3 you might not get the correct output as this system is numerically already unstable!
 *
 *  \param[in] points The points (n x 3)
 *  \param[in] L The Laplacian (n x n)
 *  \param[in] M The mass matrix (n x n)
 *  \param[in] k The smoothness degree
 */
Eigen::MatrixXd fairing(const Eigen::MatrixXd &points, const Eigen::MatrixXi &faces, const Eigen::SparseMatrix<double> &L, const int k)
{
    Eigen::MatrixXd X;
    return X;
}

/**
 * \brief The Third exercise
 *
 *  what it does:
 *  -load the pipe.obj file
 *  -construct the graph laplacian
 *  -construct the cotan laplacian
 *  -construct the normalized cotan laplacian
 *
 *  -test the fairing function for each laplacian for the values k=1, 2, 3
 *  -visualize the output
 *
 *  \param[in] filename2 The filename of the mesh
 */
void exercise3(const std::string &filename2)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    bool load_success_target = igl::readOBJ(filename2, V, F);
    if (!load_success_target)
    {
        std::cerr << "could not load file: " << filename2 << std::endl;
        return;
    }

    double min_edge_length = find_min_edge_length(V, F);
    double s = 1.0 / min_edge_length;
    scale_to_unit_cube(V, s);

    igl::viewer::Viewer viewer;

    Eigen::SparseMatrix<double> M_graph, S_graph;
    Eigen::SparseMatrix<double> L_graph = build_graph_laplacian(F, M_graph, S_graph);

    Eigen::SparseMatrix<double> M_cot, S_cot;
    Eigen::SparseMatrix<double> L_cot = build_cotan_laplacian(V, F, M_cot, S_cot, false);

    Eigen::SparseMatrix<double> M_norm, S_norm;
    Eigen::SparseMatrix<double> L_norm = build_cotan_laplacian(V, F, M_norm, S_norm, true);

    std::cout << "graph laplace, fairing, k=1" << std::endl;
    Eigen::MatrixXd X1 = fairing(V, F, L_graph, 1);
    viewer.data.set_mesh(X1, F);
    viewer.launch();

    viewer.data.clear();
    std::cout << "graph laplace, fairing, k=2" << std::endl;
    Eigen::MatrixXd X2 = fairing(V, F, L_graph, 2);
    viewer.data.set_mesh(X2, F);
    viewer.launch();


    viewer.data.clear();
    std::cout << "cot laplace, fairing, k=1" << std::endl;
    Eigen::MatrixXd X4 = fairing(V, F, L_cot, 1);
    viewer.data.set_mesh(X4, F);
    viewer.launch();

    viewer.data.clear();
    std::cout << "cot laplace, fairing, k=2" << std::endl;
    Eigen::MatrixXd X5 = fairing(V, F, L_cot, 2);
    viewer.data.set_mesh(X5, F);
    viewer.launch();


    viewer.data.clear();
    std::cout << "norm laplace, fairing, k=1" << std::endl;
    Eigen::MatrixXd X7 = fairing(V, F, L_norm, 1);
    viewer.data.set_mesh(X7, F);
    viewer.launch();

    viewer.data.clear();
    std::cout << "norm laplace, fairing, k=2" << std::endl;
    Eigen::MatrixXd X8 = fairing(V, F, L_norm, 2);
    viewer.data.set_mesh(X8, F);
    viewer.launch();

}


/**
* \brief The main function called when running this program
*
* what it does:
*  -check provided filenames
*  -run both exercises in a row
*
*  \param[in] argc The number of arguments to the binary
*  \param[in] argv The array of arguments to the binary
*/
int main(int argc, char *argv[])
{
    std::string filename1, filename2;

    if (argc == 3)
    {
        filename1 = argv[1];
        filename2 = argv[2];
    }
    else
    {
        std::cerr << "please call assignmentsheet6 like this: " << std::endl;
        std::cerr << "./bin/assignmentsheet6 data/bunny.obj data/pipe.obj" << std::endl;
        return EXIT_FAILURE;
    }

    exercise1(filename1);
    exercise2(filename1);
    exercise3(filename2);
	return EXIT_SUCCESS;
}
