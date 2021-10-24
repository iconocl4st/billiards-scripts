//
// Created by thallock on 10/19/21.
//

#include "catch2/catch.hpp"

#include "algebra/parsing.h"
#include "algebra/poly_divide.h"
#include "algebra/alg/buchberger.h"
#include "algebra/alg/f4.h"
//#include "algebra/VariableNames.h"
//#include "algebra/VarietyMatrix.h"
#include "algebra/Variety.h"
#include "extracting_from_test.h"

#include <random>

using namespace algebra::poly;


//TEST_CASE("glance_stuff", "[poly_divide]") {
//	const double ax = -20;
//	const double ay = 20;
////	const double dx = 17.759757877423255;
////	const double dy = -4.7018134271419401;
//	const double px = 0;
//	const double py = -20;
//	const double r1 = 2.26 / 2;
//	const double r2 = 2.26 / 2;
//
//	const double alpha = 2 / 7.0;
//
//	VariableNames vars;
//	vars.register_var("x");
//	vars.register_var("y"); // Location of glance
//	vars.register_var("dx");
//	vars.register_var("dy"); // Destination from glance
//	vars.register_var("s1"); // Distance along aim_line line
//	vars.register_var("s2"); // Distance along cue travel line
//
//	PolyVec x = vars("x");
//	PolyVec y = vars("y");
//	PolyVec dx = vars("dx");
//	PolyVec dy = vars("dy");
//	PolyVec s1 = vars("s1");
//	PolyVec s2 = vars("s2");
//
//	PolyVec tan_dir_x = y;
//	PolyVec tan_dir_y = -x;
//
//	PolyVec tx = x + tan_dir_x;
//	PolyVec ty = y + tan_dir_y;
//
//	PolyVec aim_dir_x = x - ax;
//	PolyVec aim_dir_y = y - ay;
//
//	PolyVec aim_x = x + s1 * (aim_dir_x - x);
//	PolyVec aim_y = y + s1 * (aim_dir_y - y);
//
//	PolyVec dest_x = x + s2 * (dx - x);
//	PolyVec dest_y = y + s2 * (dy - y);
//
//	PolyVec glance_1 = (tx - x) * alpha + (aim_x - y) * (1 - alpha) - (dx - x);
//	PolyVec glance_2 = (ty - y) * alpha + (aim_y - y) * (1 - alpha) - (dy - y);
//	PolyVec orth_req_1 = (aim_x - x) * (aim_x - tx) + (aim_y - y) * (aim_y - ty);
//	PolyVec orth_req_2 = (dest_x - x) * (dest_x - px) + (dest_y - y) * (dest_y - py);
//	PolyVec radius_1 = x.pow(2) + y.pow(2) - std::pow(r1 + r2, 2);
//	PolyVec radius_2 = (dest_x - px).pow(2) + (dest_y - py).pow(2) - std::pow(r1, 2);
//
//	glance_1 = glance_1.canon();
//	glance_2 = glance_2.canon();
//	orth_req_1 = orth_req_1.canon();
//	orth_req_2 = orth_req_2.canon();
//	radius_1 = radius_1.canon();
//	radius_2 = radius_2.canon();
//
//	std::cout << "glance_1 = " << vars(glance_1) << std::endl;
//	std::cout << "glance_2 = " << vars(glance_2) << std::endl;
//	std::cout << "orth_req_1 = " << vars(orth_req_1) << std::endl;
//	std::cout << "orth_req_2 = " << vars(orth_req_2) << std::endl;
//	std::cout << "radius_1 = " << vars(radius_1) << std::endl;
//	std::cout << "radius_2 = " << vars(radius_2) << std::endl;
//
//	std::cout << std::endl;
//
//	if (0) {
//		const auto formatter = vars.create_formatter();
//		std::vector<const PolyVec *> polynomials{&glance_1, &glance_2, &orth_req_1, &orth_req_2, &radius_1,
//												 &radius_2};
//		VarietyMatrix variety{polynomials};
//		std::cout << "Variety: " << std::endl;
//		variety.write(std::cout, formatter);
//	}
//	if (1) {
//		const auto formatter = vars.create_formatter();
//		std::vector<PolyVec> polynomials{glance_1, glance_2, orth_req_1, orth_req_2, radius_1, radius_2};
//		Variety variety;
//		variety.simplify();
////		std::cout << "Variety: " << std::endl;
////		variety.write(std::cout, formatter);
//	}
//}
//


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
	for (int i = 0; i < (int) system.system.size(); i++) {
		std::cout << i << ": " << system.system[i]->evaluate(values) << std::endl;
	}

	auto simplified = simplify::simplify_variety(system.system);

	std::cout << "After simplifying:" << std::endl;
	for (int i = 0; i < (int) simplified.size(); i++) {
		std::cout << i << ": " << simplified[i] << std::endl;
	}
	simplify::print_present_vars(simplified);

	std::cout << "Values after simplifying:" << std::endl;
	for (int i = 0; i < (int) simplified.size(); i++) {
		std::cout << i << ": " << simplified[i]->evaluate(values) << std::endl;
	}

//	std::vector<PolyPtr> basis{simplified};
#if 0
	std::vector<PolyPtr> basis{system.system};
	buchberger(basis);
#endif
#if 1
	auto basis = f4(simplified, [&](const PolyPtr& poly) {
		double value = poly->evaluate(values);
		if (std::abs(value) > 1) {
			std::cout << "Encountered polynomial " << poly << std::endl;
			throw std::runtime_error{"This polynomial is not in the ideal"};
		}
	});
//	buchberger(basis);
#endif

	std::cout << "Computed basis:" << std::endl;
	for (int i = 0; i < (int) basis.size(); i++) {
		std::cout << i << ": " << basis[i] << std::endl;
	}

	std::cout << "Basis values:" << std::endl;
	for (int i = 0; i < (int) basis.size(); i++) {
		std::cout << i << ": " << basis[i]->evaluate(values) << std::endl;
	}

	simplify::print_present_vars(basis);



//
//
//	system.variety = simplify::simplify_variety(system.variety);
//
//	std::cout << "System: " << std::endl;
//	std::cout << system << std::endl;
//
//	std::cout << "PolyVec values:" << std::endl;
//	for (int i = 0; i < (int) system.variety.polynomials.size(); i++) {
//		std::cout << i << ": " << system.variety.polynomials[i].evaluate(values) << std::endl;
//	}
}