//
// Created by thallock on 10/11/21.
//

#include "catch2/catch.hpp"

#include "algebra/parsing.h"
#include "algebra/poly_divide.h"
#include "algebra/alg/buchberger.h"
//#include "algebra/VariableNames.h"
//#include "algebra/VarietyMatrix.h"
//#include "algebra/Variety.h"
//#include "extracting_from_test.h"

#include <random>

using namespace algebra::poly;

TEST_CASE("parse poly", "[parse_polynomial]") {
	auto ideal = std::make_shared<Ideal>(cmp::Lexical);
	ideal->register_var("x[0]");
	ideal->register_var("x[1]");
	auto impl = std::make_shared<VectorIndexImpl>(ideal);
	PolyPtr p = parsing::parse_polynomial(impl, "-3x[0]^2 + 2 - 5 * x[1]");
	std::stringstream ss;
	ss << p;
	REQUIRE(ss.str() == "2 + -5 * x[1] + -3 * x[0]^2");
}

TEST_CASE("divide polynomials", "[poly_divide]") {
	auto ideal = std::make_shared<Ideal>(cmp::Lexical);
	ideal->register_var("x[0]");
	ideal->register_var("x[1]");
	auto vec_impl = std::make_shared<VectorIndexImpl>(ideal);
	auto impl = std::dynamic_pointer_cast<IndexImpl>(vec_impl);

	PolyPtr p1 = parsing::parse_polynomial(impl, "x[0] * x[1] - 1");
	PolyPtr p2 = parsing::parse_polynomial(impl, "x[1]^2 - 1");
	PolyPtr p3 = parsing::parse_polynomial(impl, "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2");
	std::cout << "p1: " << p1 << std::endl;
	std::cout << "p2: " << p2 << std::endl;
	std::cout << "p3: " << p3 << std::endl;
	Division d{p3};

	d.add_divisor(p1);
	d.add_divisor(p2);
	divide(d);

	std::cout << "Quotients:" << std::endl;
	for (const auto& quotient : d.quotients) {
		std::cout << quotient << std::endl;
	}

	std::cout << "Remainder:" << std::endl;
	std::cout << d.remainder << std::endl;

	REQUIRE(d.check());
}

TEST_CASE("divide polynomials 3", "[poly_divide]") {
	auto ideal = std::make_shared<Ideal>(cmp::Lexical);
	ideal->register_var("x[0]");
	ideal->register_var("x[1]");
	auto vec_impl = std::make_shared<VectorIndexImpl>(ideal);
	auto impl = std::dynamic_pointer_cast<IndexImpl>(vec_impl);

	PolyPtr p1 = parsing::parse_polynomial(impl, "x[1]^2 - 1");
	PolyPtr p2 = parsing::parse_polynomial(impl, "x[0] * x[1] - 1");
	PolyPtr p3 = parsing::parse_polynomial(impl, "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2");
	std::cout << "p1: " << p1 << std::endl;
	std::cout << "p2: " << p2 << std::endl;
	std::cout << "p3: " << p3 << std::endl;
	Division d{p3};
	d.add_divisor(p1);
	d.add_divisor(p2);

	divide(d);

	std::cout << "Quotients:" << std::endl;
	for (const auto& quotient : d.quotients) {
		std::cout << quotient << std::endl;
	}

	std::cout << "Remainder:" << std::endl;
	std::cout << d.remainder << std::endl;

	REQUIRE(d.check());
}


TEST_CASE("test evaluate", "[poly]") {
	auto ideal = std::make_shared<Ideal>(cmp::Lexical);
	ideal->register_var("x[0]");
	ideal->register_var("x[1]");
	auto vec_impl = std::make_shared<VectorIndexImpl>(ideal);
	auto impl = std::dynamic_pointer_cast<IndexImpl>(vec_impl);

	const auto p = parsing::parse_polynomial(impl,  "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2 - 1");
	std::cout << p << std::endl;
	std::vector<double> x{-2, 3};
	const double px = p->evaluate(x);
	REQUIRE(std::abs(px - (4 * 3 - 2 * 9 + 9 - 1)) < POLY_TOL);
}

TEST_CASE("test multiply", "[poly]") {
	auto ideal = std::make_shared<Ideal>(cmp::Lexical);
	ideal->register_var("x[0]");
	ideal->register_var("x[1]");
	auto vec_impl = std::make_shared<VectorIndexImpl>(ideal);
	auto impl = std::dynamic_pointer_cast<IndexImpl>(vec_impl);

	PolyPtr p1 = parsing::parse_polynomial(impl, "x[1]^2 - 1");
	PolyPtr p2 = parsing::parse_polynomial(impl, "x[0] * x[1] - 1");
	PolyPtr p3 = p1 * p2;

	std::vector<double> x{-2, 3};
	REQUIRE(std::abs(p1->evaluate(x) * p2->evaluate(x) - p3->evaluate(x)) < POLY_TOL);
}


TEST_CASE("test substitute", "[substitute]") {
    auto ideal = std::make_shared<Ideal>(cmp::GradedLexical);
    ideal->register_var("x[0]");
    ideal->register_var("x[1]");
    auto vec_impl = std::make_shared<VectorIndexImpl>(ideal);
    auto impl = std::dynamic_pointer_cast<IndexImpl>(vec_impl);
    PolyPtr p1 = parsing::parse_polynomial(impl, "x[0]^3 - 2 * x[0] * x[1]");
    PolyPtr p2 = parsing::parse_polynomial(impl, "x[0]^2 * x[1] - 2 * x[1]^2 + x[0]");

    simplify::
}


TEST_CASE("test buckberger", "[poly]") {
    auto ideal = std::make_shared<Ideal>(cmp::GradedLexical);
    ideal->register_var("x[0]");
    ideal->register_var("x[1]");
    auto vec_impl = std::make_shared<VectorIndexImpl>(ideal);
    auto impl = std::dynamic_pointer_cast<IndexImpl>(vec_impl);
    PolyPtr p1 = parsing::parse_polynomial(impl, "x[0]^3 - 2 * x[0] * x[1]");
    PolyPtr p2 = parsing::parse_polynomial(impl, "x[0]^2 * x[1] - 2 * x[1]^2 + x[0]");
//	PolyPtr p3 = parsing::parse_polynomial(impl, "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2");

	std::vector<PolyPtr> basis{p1, p2};

	std::cout << "Initial: " << std::endl;
	for (const auto& poly : basis) {
		std::cout << poly << std::endl;
	}
	buchberger(basis);

	std::cout << "basis: " << std::endl;
	for (const auto& poly : basis) {
		std::cout << poly << std::endl;
	}
}


// test substitute
