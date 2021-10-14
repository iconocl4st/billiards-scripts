//
// Created by thallock on 10/13/21.
//

//
// Created by thallock on 10/11/21.
//

#include "catch2/catch.hpp"

#include "algebra/f4.h"
#include "algebra/GmpPoly.h"

#include <random>
#include <iostream>

using namespace algebra::poly;

TEST_CASE("test gmp", "[gmp_arith]") {
	Ideal ideal;
	ideal.register_var("x[0]");
	ideal.register_var("x[1]");
	ideal.register_var("x[2]");

	std::string s = "-3x[0]^2 + 2 - 5 * x[1]";
	PolyDict p = parse::parse_polynomial(ideal, s);
	std::stringstream ss;
	ss << p;
	REQUIRE(ss.str() == "2 + -3 * (2) x[0]^2 + -5 * (256) x[1]");
}

TEST_CASE("gmp test evaluate", "[poly]") {
	Ideal ideal;
	ideal.register_var("x[0]");
	ideal.register_var("x[1]");
	std::string s = "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2 - 1";
	const PolyDict p = parse::parse_polynomial(ideal, s);
	std::vector<double> x{-2, 3};
	const double px = p.evaluate(x);
	REQUIRE(std::abs(px - (4 * 3 - 2 * 9 + 9 - 1)) < YET_ANOTHER_TOL);
}

TEST_CASE("gmp test multiply", "[poly]") {
	Ideal ideal;
	ideal.register_var("x[0]");
	ideal.register_var("x[1]");
	std::string s1 = "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2 + 2";
	std::string s2 = "x[0]^3 * x[1] + x[0] - 3";
	const PolyDict p1 = parse::parse_polynomial(ideal, s1);
	const PolyDict p2 = parse::parse_polynomial(ideal, s2);
	const PolyDict p3 = p1 * p2;
	std::vector<double> x{-2, 3};
	REQUIRE(std::abs(p3.evaluate(x) - p1.evaluate(x) * p2.evaluate(x)) < YET_ANOTHER_TOL);
}

TEST_CASE("f4", "[poly]") {
	Ideal ideal;
	ideal.register_var("x");
	ideal.register_var("y");
	ideal.register_var("z");
	std::string s1 = "x[0] * x[2] - x[1]^2";
	std::string s2 = "x[0]^3 - x[2]^2";
	// page 97
	const PolyDict p1 = parse::parse_polynomial(ideal, s1);
	const PolyDict p2 = parse::parse_polynomial(ideal, s2);

	std::vector<PolyDict> polys{p1, p2};

	std::cout << "Polynomial generators:" << std::endl;
	for (const auto& b : polys) {
		std::cout << "\t" << b << std::endl;
	}

	f4(polys);

	std::cout << "Computed basis:" << std::endl;
	for (const auto& b : polys) {
		std::cout << "\t" << b << std::endl;
	}
}
