#include <iostream>
#include <vector>
#include <RcppEigen.h>
#include "include/mvn_stats.h"
#include "include/mvn_stats_interpoint.h"
#include "include/mvn_stats_approx.h"
#include "include/kronecker.h"

using namespace mvn;
using namespace matching;

// Helper function to convert clustering from R to Clustering object
Clustering convert_clustering_from_r(const std::vector<int>& clustering) {
    return Clustering(clustering);
}

// Helper function to create mahalanobis_distances
std::pair<mahalanobis_distances, mahalanobis_distances> create_mahalanobis_distances(
    const Eigen::MatrixXd& matrix, const Eigen::MatrixXd& S, const Eigen::VectorXd& mean) {
    auto matrix_ptr = std::make_shared<Eigen::MatrixXd>(matrix);
    mahalanobis_distances distances_interpoint(matrix_ptr, S, mean);
    mahalanobis_distances distances_approx(matrix_ptr, S, mean);
    return {distances_interpoint, distances_approx};
}

// Helper function to convert Clustering to std::vector<std::vector<int>>
std::vector<std::vector<int>> convert_clustering(const Clustering& clst) {
    std::vector<std::vector<int>> clusters(clst.size());
    for (size_t i = 0; i < clst.size(); ++i) {
        clusters[i] = clst.elements(i);
    }
    return clusters;
}

void preprocess(const Eigen::MatrixXd& matrix, const Clustering& clst, const std::string& filename) {
    kronecker_preprocessor preprocessor(std::make_shared<Eigen::MatrixXd>(matrix), convert_clustering(clst), filename);
    preprocessor.process(1, 1, 4);
}

// Helper function to process stats
Rcpp::List process_stats(const Eigen::MatrixXd& matrix, const Eigen::MatrixXd& S, const Eigen::VectorXd& mean, const Clustering& clst) {
    auto matrix_transposed = matrix.transpose();
    auto [distances_interpoint, distances_approx] = create_mahalanobis_distances(matrix_transposed, S, mean);
    distances_interpoint.calculate_interpoint();

    for (size_t i = 0; i < matrix.rows() ; ++i) {
        for (size_t j = 0; j < matrix.rows(); ++j) {
            Rcpp::Rcout << distances_interpoint.interpoint_distance(i, j) << " ";

        }
        Rcpp::Rcout << std::endl;
    }

    auto clusters = convert_clustering(clst);

    std::string filename = "kronecker_data.bin";
    preprocess(matrix, clst, filename);

    double beta = 0.3;

    mvn_stats_interpoint interpoint_stats(clst.size());
    interpoint_stats.init(distances_interpoint, clst, beta);
    interpoint_stats.init_pairwise(distances_interpoint, clst, beta);

    mvn_stats_approx approx_stats(filename, clst.size());
    approx_stats.init(distances_approx, clst, beta);
    approx_stats.init_pairwise(distances_approx, clst, beta);

    std::vector<double> interpoint_centered, approx_centered;
    for (size_t i = 0; i < clst.size(); ++i) {
        interpoint_centered.push_back(interpoint_stats.centered_stat(i));
        approx_centered.push_back(approx_stats.centered_stat(i));
    }

    std::vector<std::vector<double>> interpoint_pairwise(clst.size()), approx_pairwise(clst.size());
    for (size_t i = 0; i < clst.size(); ++i) {
        for (size_t j = 0; j < clst.size(); ++j) {
            interpoint_pairwise[i].push_back(interpoint_stats.pairwise_stat(i, j));
            approx_pairwise[i].push_back(approx_stats.pairwise_stat(i, j));
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("interpoint_centered") = interpoint_centered,
        Rcpp::Named("approx_centered") = approx_centered,
        Rcpp::Named("interpoint_pairwise") = interpoint_pairwise,
        Rcpp::Named("approx_pairwise") = approx_pairwise
    );
}

// Combined function to handle both centered and pairwise stats
Rcpp::List run_mvn_stats_tests_combined(const Eigen::MatrixXd& S, const Eigen::VectorXd& mean, const Eigen::MatrixXd& matrix, const std::vector<int>& clustering) {
    Clustering clst = convert_clustering_from_r(clustering);
    // print test info using Rcpp
    Rcpp::Rcout << "Running mvn stats tests with " << clst.size() << " clusters." << std::endl;
    Rcpp::Rcout << "Matrix size: " << matrix.rows() << " x " << matrix.cols() << std::endl;
    Rcpp::Rcout.flush();
    return process_stats(matrix, S, mean, clst);
}

// [[Rcpp::export]]
Rcpp::List rcpp_run_mvn_stats_tests_combined(const Eigen::MatrixXd& S, const Eigen::VectorXd& mean, const Eigen::MatrixXd& matrix, const std::vector<int>& clustering) {
    return run_mvn_stats_tests_combined(S, mean, matrix, clustering);
}
