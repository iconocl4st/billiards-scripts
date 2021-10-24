//
// Created by thallock on 10/21/21.
//


#include "catch2/catch.hpp"

#include <random>
#include <iostream>

#include <Eigen/SparseCore>
#include <Eigen/Sparse>

#include "algebra/alg/sparse_jordan.h"


TEST_CASE("test sparse reduction", "[reduce]") {
	int nrows = 6;
	int ncols = 10;

	std::vector<Eigen::Triplet<double>> triples;
	triples.reserve(nrows * ncols);

	Eigen::SparseMatrix<double, Eigen::RowMajor> M{nrows, ncols};

	std::random_device rd;
	std::default_random_engine e1(rd());
	std::uniform_real_distribution<double> c_dist(-10, 10);
	std::uniform_real_distribution<double> b_dist(0, 1);

	for (int r = 0; r < nrows; r++) {
		for (int c = 0; c < ncols; c++) {
			if (b_dist(e1) < 0.9) {
				continue;
			}
			triples.emplace_back(r, c, c_dist(e1));
		}
	}
	M.setFromTriplets(triples.begin(), triples.end());

	std::cout << "Matrix:" << std::endl;
	std::cout << M.toDense() << std::endl;

	algebra::poly::sparse_jordan(M);

	std::cout << "Reduced matrix:" << std::endl;
	std::cout << M.toDense() << std::endl;
}