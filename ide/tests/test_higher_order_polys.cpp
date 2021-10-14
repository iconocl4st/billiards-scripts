//
// Created by thallock on 10/5/21.
//
//
// Created by thallock on 10/3/21.
//

#include "catch2/catch.hpp"

#include "math/high_order_polys.h"


/*
import itertools


def print_coeffs(num_roots):
	roots = ['r' + str(i) for i in range(1, num_roots + 1)]
	leading = 'c' + str(num_roots)
	for ord in range(1, num_roots + 1):
		print('\tconst double c' + str(num_roots - ord) + ' = ' + ('+' if ord % 2 == 0 else '-') + leading + ' * (' + ' + '.join(
			'*'.join(factors)
			for factors in itertools.combinations(roots, r=ord)) + ');')

print_coeffs(6)


def print_poly(degree):
	print('\t const auto eval = [&](const double x) { return ' + ' + '.join(
		'c' + str(idx) + ' * std::pow(x, ' + str(idx) + ')'
		if idx > 1 else (
			'c1 * x' if idx == 1 else 'c0'
		)
		for idx in range(degree + 1)
	) + '; };')

print_poly(6)

 */

TEST_CASE("Solve higher order polynomials happy case", "[compute_roots]") {
	const double r1 = 3;
	const double r2 = 4;
	const double r3 = -2;
	const double r4 = 17;
	const double r5 = M_PI;
	const double r6 = 2;

	const double c6 = 2;
	const double c5 = -c6 * (r1 + r2 + r3 + r4 + r5 + r6);
	const double c4 = +c6 * (r1*r2 + r1*r3 + r1*r4 + r1*r5 + r1*r6 + r2*r3 + r2*r4 + r2*r5 + r2*r6 + r3*r4 + r3*r5 + r3*r6 + r4*r5 + r4*r6 + r5*r6);
	const double c3 = -c6 * (r1*r2*r3 + r1*r2*r4 + r1*r2*r5 + r1*r2*r6 + r1*r3*r4 + r1*r3*r5 + r1*r3*r6 + r1*r4*r5 + r1*r4*r6 + r1*r5*r6 + r2*r3*r4 + r2*r3*r5 + r2*r3*r6 + r2*r4*r5 + r2*r4*r6 + r2*r5*r6 + r3*r4*r5 + r3*r4*r6 + r3*r5*r6 + r4*r5*r6);
	const double c2 = +c6 * (r1*r2*r3*r4 + r1*r2*r3*r5 + r1*r2*r3*r6 + r1*r2*r4*r5 + r1*r2*r4*r6 + r1*r2*r5*r6 + r1*r3*r4*r5 + r1*r3*r4*r6 + r1*r3*r5*r6 + r1*r4*r5*r6 + r2*r3*r4*r5 + r2*r3*r4*r6 + r2*r3*r5*r6 + r2*r4*r5*r6 + r3*r4*r5*r6);
	const double c1 = -c6 * (r1*r2*r3*r4*r5 + r1*r2*r3*r4*r6 + r1*r2*r3*r5*r6 + r1*r2*r4*r5*r6 + r1*r3*r4*r5*r6 + r2*r3*r4*r5*r6);
	const double c0 = +c6 * (r1*r2*r3*r4*r5*r6);

	const auto eval = [&](const double x) { return c0 + c1 * x + c2 * std::pow(x, 2) + c3 * std::pow(x, 3) + c4 * std::pow(x, 4) + c5 * std::pow(x, 5) + c6 * std::pow(x, 6); };

	std::vector<double> coefficients{c0, c1, c2, c3, c4, c5, c6};
	std::list<double> roots;
	billiards::shots::math::compute_roots(coefficients, roots);
	for (const double root : roots) {
		REQUIRE(std::abs(eval(root)) < 1e-4);
	}
	REQUIRE(roots.size() == 6);
}


TEST_CASE("Solve higher order polynomials repeated roots", "[compute_roots]") {
	const double r1 = 3;
	const double r2 = 3;
	const double r3 = 3;
	const double r4 = 17;
	const double r5 = M_PI;
	const double r6 = 2;

	const double c6 = 2;
	const double c5 = -c6 * (r1 + r2 + r3 + r4 + r5 + r6);
	const double c4 = +c6 * (r1*r2 + r1*r3 + r1*r4 + r1*r5 + r1*r6 + r2*r3 + r2*r4 + r2*r5 + r2*r6 + r3*r4 + r3*r5 + r3*r6 + r4*r5 + r4*r6 + r5*r6);
	const double c3 = -c6 * (r1*r2*r3 + r1*r2*r4 + r1*r2*r5 + r1*r2*r6 + r1*r3*r4 + r1*r3*r5 + r1*r3*r6 + r1*r4*r5 + r1*r4*r6 + r1*r5*r6 + r2*r3*r4 + r2*r3*r5 + r2*r3*r6 + r2*r4*r5 + r2*r4*r6 + r2*r5*r6 + r3*r4*r5 + r3*r4*r6 + r3*r5*r6 + r4*r5*r6);
	const double c2 = +c6 * (r1*r2*r3*r4 + r1*r2*r3*r5 + r1*r2*r3*r6 + r1*r2*r4*r5 + r1*r2*r4*r6 + r1*r2*r5*r6 + r1*r3*r4*r5 + r1*r3*r4*r6 + r1*r3*r5*r6 + r1*r4*r5*r6 + r2*r3*r4*r5 + r2*r3*r4*r6 + r2*r3*r5*r6 + r2*r4*r5*r6 + r3*r4*r5*r6);
	const double c1 = -c6 * (r1*r2*r3*r4*r5 + r1*r2*r3*r4*r6 + r1*r2*r3*r5*r6 + r1*r2*r4*r5*r6 + r1*r3*r4*r5*r6 + r2*r3*r4*r5*r6);
	const double c0 = +c6 * (r1*r2*r3*r4*r5*r6);

	const auto eval = [&](const double x) { return c0 + c1 * x + c2 * std::pow(x, 2) + c3 * std::pow(x, 3) + c4 * std::pow(x, 4) + c5 * std::pow(x, 5) + c6 * std::pow(x, 6); };

	std::vector<double> coefficients{c0, c1, c2, c3, c4, c5, c6};
	std::list<double> roots;
	billiards::shots::math::compute_roots(coefficients, roots);
	for (const double root : roots) {
		REQUIRE(std::abs(eval(root)) < 1e-4);
	}
	REQUIRE(roots.size() == 4);
}

TEST_CASE("Get roots of linear", "[compute_roots]") {
	std::vector<double> coefficients{3, 2};

	std::list<double> roots;
	billiards::shots::math::compute_roots(coefficients, roots);
	for (const double root : roots) {
		REQUIRE(std::abs(root + 3 / 2.0) < 1e-4);
	}
	REQUIRE(roots.size() == 1);
}

TEST_CASE("Get roots of quadratic", "[compute_roots]") {
	std::vector<double> coefficients{6, 5, 1};

	std::list<double> roots;
	billiards::shots::math::compute_roots(coefficients, roots);
	for (const double root : roots) {
		REQUIRE((
			std::abs(root + 2) < 1e-4 ||
			std::abs(root + 3) < 1e-4
		));
	}
	REQUIRE(roots.size() == 2);
}

TEST_CASE("Solve all zeros", "[compute_roots]") {
	std::vector<double> coefficients{0, 0, 0, 0, 0};
	std::list<double> roots;
	billiards::shots::math::compute_roots(coefficients, roots);
	for (const double root : roots) {
		REQUIRE(std::abs(root) < 1e-4);
	}
	REQUIRE(roots.size() == 1);
}

TEST_CASE("Solve higher order polynomials impossible", "[compute_roots]") {
	std::vector<double> coefficients{1, 0, 0, 0, 0};
	std::list<double> roots;
	billiards::shots::math::compute_roots(coefficients, roots);
	REQUIRE(roots.empty());
}

TEST_CASE("Solve higher order polynomials imaginary roots", "[compute_roots]") {
	std::vector<double> coefficients{1, 0, 1};
	std::list<double> roots;
	billiards::shots::math::compute_roots(coefficients, roots);
	REQUIRE(roots.empty());
}
