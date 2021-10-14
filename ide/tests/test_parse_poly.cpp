//
// Created by thallock on 10/11/21.
//

#include "catch2/catch.hpp"

#include "algebra/parsing.h"
#include "algebra/buchberger.h"
#include "algebra/poly_divide.h"
#include "algebra/VariableNames.h"
#include "algebra/VarietyMatrix.h"
#include "algebra/Variety.h"
#include "extracting_from_test.h"

#include <random>

using namespace algebra::poly;

TEST_CASE("parse poly", "[parse_polynomial]") {
	std::string s = "-3x[0]^2 + 2 - 5 * x[1]";
	Polynomial p = parsing::parse_polynomial(s, 3);
	std::stringstream ss;
	ss << p;
	REQUIRE(ss.str() == "-3 * x[0]^2  + 2 + -5");
}

TEST_CASE("collect polynomial", "[poly collect]") {
	std::string s = "-3x[0]^2 + 2 + 5 * x[0]^2";
	Polynomial p1 = parsing::parse_polynomial(s, 3);
	Polynomial p2 = p1.collect();
	std::cout << p2 << std::endl;
	std::stringstream ss;
	ss << p2;
	REQUIRE(ss.str() == "2 * x[0]^2  + 2");
}

TEST_CASE("sort polynomial", "[poly collect]") {
	std::string s = "-3x[1]^2 + 2 + 3 * x[0] * x[2] + 5 * x[0]^2";
	Polynomial p = parsing::parse_polynomial(s, 3);
	p.sort();
	std::stringstream ss;
	ss << p;
	REQUIRE(ss.str() == "5 * x[0]^2  + 3 * x[0] x[2]  + -3 * x[1]^2  + 2");
}


TEST_CASE("divide polynomials", "[poly_divide]") {
	std::string s1 = "x[0] * x[1] - 1";
	std::string s2 = "x[1]^2 - 1";
	std::string s3 = "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2";
	const Polynomial p1 = parsing::parse_polynomial(s1, 2);
	const Polynomial p2 = parsing::parse_polynomial(s2, 2);
	const Polynomial p3 = parsing::parse_polynomial(s3, 2);
	std::cout << "p1: "  << p1 << std::endl;
	std::cout << "p2: "  << p2 << std::endl;
	std::cout << "p3: " << p3 << std::endl;
	std::vector<Polynomial> divisors{p1, p2};
	Division d{divisors, p3};
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
	std::string s1 = "x[1]^2 - 1";
	std::string s2 = "x[0] * x[1] - 1";
	std::string s3 = "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2";
	const Polynomial p1 = parsing::parse_polynomial(s1, 2);
	const Polynomial p2 = parsing::parse_polynomial(s2, 2);
	const Polynomial p3 = parsing::parse_polynomial(s3, 2);
	std::cout << "p1: "  << p1 << std::endl;
	std::cout << "p2: "  << p2 << std::endl;
	std::cout << "p3: " << p3 << std::endl;
	std::vector<Polynomial> divisors{p1, p2};
	Division d{divisors, p3};
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
	std::string s = "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2 - 1";
	const Polynomial p = parsing::parse_polynomial(s, 2);
	std::vector<double> x{-2, 3};
	const double px = p.evaluate(x);
	REQUIRE(std::abs(px - (4 * 3 - 2 * 9 + 9 - 1)) < POLY_TOL);
}

TEST_CASE("test multiply", "[poly]") {
	std::string s1 = "x[0]^2 * x[1] + x[0] * x[1]^2 + x[1]^2 + 2";
	std::string s2 = "x[0]^3 * x[1] + x[0] - 3";
	const Polynomial p1 = parsing::parse_polynomial(s1, 2);
	const Polynomial p2 = parsing::parse_polynomial(s2, 2);
	const Polynomial p3 = p1 * p2;
	std::vector<double> x{-2, 3};
	REQUIRE(std::abs(p3.evaluate(x) - p1.evaluate(x) * p2.evaluate(x)) < POLY_TOL);
}




TEST_CASE("glance_stuff", "[poly_divide]") {
	const double ax = -20;
	const double ay = 20;
//	const double dx = 17.759757877423255;
//	const double dy = -4.7018134271419401;
	const double px = 0;
	const double py = -20;
	const double r1 = 2.26 / 2;
	const double r2 = 2.26 / 2;

	const double alpha = 2 / 7.0;

	VariableNames vars;
	vars.register_var("x");
	vars.register_var("y"); // Location of glance
	vars.register_var("dx");
	vars.register_var("dy"); // Destination from glance
	vars.register_var("s1"); // Distance along aim_line line
	vars.register_var("s2"); // Distance along cue travel line

	Polynomial x = vars("x");
	Polynomial y = vars("y");
	Polynomial dx = vars("dx");
	Polynomial dy = vars("dy");
	Polynomial s1 = vars("s1");
	Polynomial s2 = vars("s2");

	Polynomial tan_dir_x = y;
	Polynomial tan_dir_y = -x;

	Polynomial tx = x + tan_dir_x;
	Polynomial ty = y + tan_dir_y;

	Polynomial aim_dir_x = x - ax;
	Polynomial aim_dir_y = y - ay;

	Polynomial aim_x = x + s1 * (aim_dir_x - x);
	Polynomial aim_y = y + s1 * (aim_dir_y - y);

	Polynomial dest_x = x + s2 * (dx - x);
	Polynomial dest_y = y + s2 * (dy - y);

	Polynomial glance_1 = (tx - x) * alpha + (aim_x - y) * (1 - alpha) - (dx - x);
	Polynomial glance_2 = (ty - y) * alpha + (aim_y - y) * (1 - alpha) - (dy - y);
	Polynomial orth_req_1 = (aim_x - x) * (aim_x - tx) + (aim_y - y) * (aim_y - ty);
	Polynomial orth_req_2 = (dest_x - x) * (dest_x - px) + (dest_y - y) * (dest_y - py);
	Polynomial radius_1 = x.pow(2) + y.pow(2) - std::pow(r1 + r2, 2);
	Polynomial radius_2 = (dest_x - px).pow(2) + (dest_y - py).pow(2) - std::pow(r1, 2);

	glance_1 = glance_1.canon();
	glance_2 = glance_2.canon();
	orth_req_1 = orth_req_1.canon();
	orth_req_2 = orth_req_2.canon();
	radius_1 = radius_1.canon();
	radius_2 = radius_2.canon();

	std::cout << "glance_1 = " << vars(glance_1) << std::endl;
	std::cout << "glance_2 = " << vars(glance_2) << std::endl;
	std::cout << "orth_req_1 = " << vars(orth_req_1) << std::endl;
	std::cout << "orth_req_2 = " << vars(orth_req_2) << std::endl;
	std::cout << "radius_1 = " << vars(radius_1) << std::endl;
	std::cout << "radius_2 = " << vars(radius_2) << std::endl;

	std::cout << std::endl;

	if (0) {
		const auto formatter = vars.create_formatter();
		std::vector<const Polynomial *> polynomials{&glance_1, &glance_2, &orth_req_1, &orth_req_2, &radius_1,
													&radius_2};
		VarietyMatrix variety{polynomials};
		std::cout << "Variety: " << std::endl;
		variety.write(std::cout, formatter);
	}
	if (1) {
		const auto formatter = vars.create_formatter();
		std::vector<Polynomial> polynomials{glance_1, glance_2, orth_req_1, orth_req_2, radius_1, radius_2};
		Variety variety;
		variety.simplify();
//		std::cout << "Variety: " << std::endl;
//		variety.write(std::cout, formatter);
	}
}

TEST_CASE("create solution", "[poly_divide]") {
	std::random_device r;
	std::default_random_engine e1(r());
	std::uniform_real_distribution<double> ball_location_dist(-10, 10);
	std::uniform_real_distribution<double> radius_dist(0, 2);
	std::uniform_real_distribution<double> zo_dist(-1, 1);

	double ax = ball_location_dist(e1);
	double ay = ball_location_dist(e1);
	double ox = ball_location_dist(e1);
	double oy = ball_location_dist(e1);
	double r1 = radius_dist(e1);
	double r2 = radius_dist(e1);
	double gx = ox + (r1 + r2) * zo_dist(e1);

	double y_offset = std::sqrt(std::pow(r1 + r2, 2) - std::pow(gx - ox, 2));
	double gy = oy + y_offset;
	if ((gx - ox) * ax + gy * ay < 0) {
		gy = oy - y_offset;
	}
	const double aim_dir_x = gx - ax;
	const double aim_dir_y = gy - ay;
	double tan_dir_x = (gy - oy);
	double tan_dir_y = -(gx - ox);
	if (aim_dir_x * tan_dir_x + aim_dir_y * tan_dir_y < 0 && false) {
		tan_dir_x = -(gy - oy);
		tan_dir_y = (gx - ox);
	}
	double tx = gx + tan_dir_x;
	double ty = gy + tan_dir_y;
	// Want:
	// (x + s * aim_dir_x - x) * (x + s * aim_dir_x - tx) + (y + s * aim_dir_y - y) * (y + s * aim_dir_y - ty) == 0
	// (s * aim_dir_x) * (s * aim_dir_x - tan_dir_x) + (s * aim_dir_y) * (s * aim_dir_y - tan_dir_y) == 0
	// aim_dir_x * (s * aim_dir_x - tan_dir_x) + aim_dir_y * (s * aim_dir_y - tan_dir_y) == 0
	// s * aim_dir_x * aim_dir_x - aim_dir_x * tan_dir_x + s * aim_dir_y * aim_dir_y - aim_dir_y * tan_dir_y == 0
	// s * aim_dir_x * aim_dir_x + s * aim_dir_y * aim_dir_y  ==  aim_dir_x * tan_dir_x + aim_dir_y * tan_dir_y
	// s ==  (aim_dir_x * tan_dir_x + aim_dir_y * tan_dir_y) / (aim_dir_x * aim_dir_x + aim_dir_y * aim_dir_y)
	double s1 = (aim_dir_x * tan_dir_x + aim_dir_y * tan_dir_y) / (aim_dir_x * aim_dir_x + aim_dir_y * aim_dir_y);
	double aim_x = gx + s1 * aim_dir_x;
	double aim_y = gy + s1 * aim_dir_y;
	std::cout << "Generation ensures this is 0: " << (aim_x - gx) * (aim_x - tx) + (aim_y - gy) * (aim_y - ty) << std::endl;
	double alpha = 2 / 7.0;
	double dx = tx * (1 - alpha) + aim_x * alpha;
	double dy = ty * (1 - alpha) + aim_y * alpha;
	double s2 = std::abs(ball_location_dist(e1));
//	double norm_dist = std::sqrt(std::pow(dx - gx, 2) + std::pow(dy - gy, 2));
	double dest_x = gx + s2 * (dx - gx);
	double dest_y = gy + s2 * (dy - gy);
	double pocket_dir_x;
	double pocket_dir_y;

	if (zo_dist(e1) < 0) {
		pocket_dir_x = (dy - gy);
		pocket_dir_y = -(dx - gx);
	} else {
		pocket_dir_x = -(dy - gy);
		pocket_dir_y = (dx - gx);
	}
	double dest_dir_norm = std::sqrt(pocket_dir_x * pocket_dir_x + pocket_dir_y * pocket_dir_y);
	pocket_dir_x *= r1 / dest_dir_norm;
	pocket_dir_y *= r1 / dest_dir_norm;
	double px = dest_x + pocket_dir_x;
	double py = dest_y + pocket_dir_y;

	std::cout << "glance: (" << gx << ", " << gy << ")" << std::endl;
	std::cout << "cue: (" << ax << ", " << ay << ")" << std::endl;
	std::cout << "linear comb: (" << dx << ", " << dy << ")" << std::endl;
	std::cout << "tangent: (" << tx << ", " << ty << ")" << std::endl;
	std::cout << "pocket: (" << px << ", " << py << ")" << std::endl;
	std::cout << "dest: (" << dest_x << ", " << dest_y << ")" << std::endl;
	std::cout << "aim: (" << aim_x << ", " << aim_y << ")" << std::endl;
	std::cout << "object: (" << ox << ", " << oy << ")" << std::endl;

	// Variables (10): s_13_x=x[0], s_13_y=x[1], t_13_x=x[2], t_13_y=x[3], d_13_x=x[4], d_13_y=x[5], u_13_1=x[6], u_13_2=x[7], t_14_x=x[8], t_14_y=x[9]
	std::vector<double> values{gx, gy, gx, gy, dx, dy, s1, s2, dest_x, dest_y};

	billiards::shots::math::LocationsSystem system;
	system.src_x = ax;
	system.src_y = ay;
	system.begin_index = 13;

	billiards::shots::math::register_rolling_glance_vars(system, system.begin_index);
	billiards::shots::math::add_pocket_vars(system, system.begin_index + 1);

	billiards::shots::math::add_glance(system, r1, r2, ox, oy, system.begin_index);
	billiards::shots::math::add_pocket(system, r1, px, py, system.begin_index + 1);

	std::cout << "System: " << std::endl;
	std::cout << system << std::endl;

	std::cout << "Polynomial values:" << std::endl;
	for (int i = 0; i < (int) system.variety.polynomials.size(); i++) {
		std::cout << i << ": " << system.variety.polynomials[i].evaluate(values) << std::endl;
	}

	std::vector<Polynomial> basis{system.variety.polynomials};
//	buchberger(basis);


	system.variety = simplify::simplify_variety(system.variety);

	std::cout << "System: " << std::endl;
	std::cout << system << std::endl;

	std::cout << "Polynomial values:" << std::endl;
	for (int i = 0; i < (int) system.variety.polynomials.size(); i++) {
		std::cout << i << ": " << system.variety.polynomials[i].evaluate(values) << std::endl;
	}
}