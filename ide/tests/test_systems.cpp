//
// Created by thallock on 10/4/21.
//

#include "catch2/catch.hpp"

#include "billiards_common/math/quad_systems.h"
#include "SolutionsTracker.h"

using namespace billiards::math;

///////////////////////////////////////////////////////
// solve_20
///////////////////////////////////////////////////////
TEST_CASE("system 02 happy case", "[solve_20]") {
	const double c = 2;
	const double r1 = 2;
	const double r2 = -3;
	const double a00 = c * r1 * r2;
	const double a10 = -c * (r1 + r2);
	const double a20 = c;
	const double default_x = 13;
	const double default_y = 15;
	billiards::test::SolutionsTracker2d tracker;
	solve_20(a00, a10, a20, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 2);
	REQUIRE(tracker.contains(std::pair<double, double>{r1, default_y}));
	REQUIRE(tracker.contains(std::pair<double, double>{r2, default_y}));
	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(a00 + a10 * p.first + a20 * p.first * p.first) < LARGER_TOL);
	}
}


///////////////////////////////////////////////////////
// solve_21
///////////////////////////////////////////////////////
TEST_CASE("system 21 happy case", "[solve_21]") {
	const double a00 = -12;
	const double a10 = 2;
	const double a01 = 2;
	const double a20 = 4;
	const double default_x = 13;
	const double default_y = 15;
	billiards::test::SolutionsTracker2d tracker;
	solve_21(a00, a10, a01, a20, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE((tracker.sols->at(0).first - default_x) < LARGER_TOL);
	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(a00 + a10 * p.first + a01 * p.second + a20 * p.first * p.first) < LARGER_TOL);
	}
}
TEST_CASE("system 21 single variable", "[solve_21]") {
	const double a00 = -12;
	const double a10 = 2;
	const double a01 = 0;
	const double a20 = 4;
	const double default_x = 13;
	const double default_y = 15;
	billiards::test::SolutionsTracker2d tracker;
	solve_21(a00, a10, a01, a20, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 2);
	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(a00 + a10 * p.first + a01 * p.second + a20 * p.first * p.first) < LARGER_TOL);
	}
}
TEST_CASE("system 21 no solutions", "[solve_21]") {
	const double a00 = 1;
	const double a10 = 0;
	const double a01 = 0;
	const double a20 = 0;
	const double default_x = 13;
	const double default_y = 15;
	billiards::test::SolutionsTracker2d tracker;
	solve_21(a00, a10, a01, a20, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 0);
}


///////////////////////////////////////////////////////
// solve_22
///////////////////////////////////////////////////////
TEST_CASE("system 22 happy case", "[solve_22]") {
	const double a00 = -5;
	const double a10 = 0;
	const double a01 = 0;
	const double a20 = 0;
	const double a11 = 3;
	const double a02 = 3;
	const double default_x = 13;
	const double default_y = 15;

	billiards::test::SolutionsTracker2d tracker;
	solve_22(a00, a10, a01, a20, a11, a02, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 2);
	REQUIRE((tracker.sols->at(0).first - default_x) < LARGER_TOL);
	REQUIRE((tracker.sols->at(1).first - default_x) < LARGER_TOL);
	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(
			a00 +
			a10 * p.first + a01 * p.second +
			a20 * p.first * p.first + a11 * p.first * p.second + a02 * p.second * p.second) < LARGER_TOL);
	}
}


TEST_CASE("system 22 discriminant has no roots and default_x works", "[solve_22]") {
	/*
ass = {'a01': 3, 'a20': -4, 'a11': 2, 'a02': 4, 'x': 0.13}
a10_sol = solve(disc.substitute(**ass) == 1, a10)[0].right()
a00_sol = solve(disc_disc.substitute(**ass, a10=a10_sol) == 1, a00)[0].right()
print('\tconst double a00 = ' + convert_powers(a00_sol) + ';')
print('\tconst double a10 = ' + convert_powers(a10_sol) + ';')

print('\t\tconst double disc_disc = ' + convert_powers(disc_disc) + ';')
print('\t\tconst double disc = ' + convert_powers(disc) + ';')
*/
	const double a00 = -13/ (double) 1600*sqrt(273) + 17127/ (double) 40000;
	const double a10 = -100/ (double) 13*a00 + 26773/ (double) 5200;
	const double a01 = 3;
	const double a20 = -4;
	const double a11 = 2;
	const double a02 = 4;
	const double default_x = 0.13;
	const double default_y = 0.15;

	{
		const double x = default_x;
		const double disc_disc = 16*std::pow(a02, 2)*std::pow(a10, 2) - 16*a01*a02*a10*a11 + 16*a00*a02*std::pow(a11, 2) + 16*std::pow(a01, 2)*a02*a20 - 64*a00*std::pow(a02, 2)*a20;
		const double disc = (std::pow(a11, 2) - 4*a02*a20)*std::pow(x, 2) + std::pow(a01, 2) - 4*a00*a02 - 2*(2*a02*a10 - a01*a11)*x;
		REQUIRE(std::abs(disc_disc - 1) < LARGER_TOL);
		REQUIRE(std::abs(disc - 1) < LARGER_TOL);
	}

	billiards::test::SolutionsTracker2d tracker;
	solve_22(a00, a10, a01, a20, a11, a02, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 2);
	REQUIRE((tracker.sols->at(0).first - default_x) < LARGER_TOL);
	REQUIRE((tracker.sols->at(1).first - default_x) < LARGER_TOL);
	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(
			a00 +
			a10 * p.first + a01 * p.second +
			a20 * p.first * p.first + a11 * p.first * p.second + a02 * p.second * p.second) < LARGER_TOL);
	}
}

TEST_CASE("system 22 discriminant has no roots and default_x does not work", "[solve_22]") {
	/*
ass = {'a01': 3, 'a20': 4, 'a11': 2, 'a02': 4, 'x': 0.13}
a10_sol = solve(disc.substitute(**ass) == -1, a10)[0].right()
a00_sol = solve(disc_disc.substitute(**ass, a10=a10_sol) == -1, a00)[0].right()
print('\tconst double a00 = ' + convert_powers(a00_sol) + ';')
print('\tconst double a10 = ' + convert_powers(a10_sol) + ';')
*/
	const double a00 = -13/ (double) 1600*sqrt(239) + 5507/ (double) 8000;
	const double a10 = -100/ (double) 13*a00 + 5273/ (double) 1040;
	const double a01 = 3;
	const double a20 = 4;
	const double a11 = 2;
	const double a02 = 4;
	const double default_x = 0.13;
	const double default_y = 0.15;

	{
		const double x = default_x;
		const double disc_disc = 16*std::pow(a02, 2)*std::pow(a10, 2) - 16*a01*a02*a10*a11 + 16*a00*a02*std::pow(a11, 2) + 16*std::pow(a01, 2)*a02*a20 - 64*a00*std::pow(a02, 2)*a20;
		const double disc = (std::pow(a11, 2) - 4*a02*a20)*std::pow(x, 2) + std::pow(a01, 2) - 4*a00*a02 - 2*(2*a02*a10 - a01*a11)*x;
		REQUIRE(std::abs(disc_disc + 1) < LARGER_TOL);
		REQUIRE(std::abs(disc + 1) < LARGER_TOL);
	}

	billiards::test::SolutionsTracker2d tracker;
	solve_22(a00, a10, a01, a20, a11, a02, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 0);
}

TEST_CASE("system 22 discriminant has one root and default_x does not work", "[solve_22]") {
	/*
ass = {'a01': 3, 'a20': 4, 'a11': 2, 'a02': 4, 'x': 0.13}
a10_sol = solve(disc.substitute(**ass) == -1, a10)[0].right()
a00_sol = solve(disc_disc.substitute(**ass, a10=a10_sol) == 0, a00)[0].right()
print('\tconst double a00 = ' + convert_powers(a00_sol) + ';')
print('\tconst double a10 = ' + convert_powers(a10_sol) + ';')

*/
	const double a00 = -13/ (double) 400*sqrt(15) + 5507/ (double) 8000;
	const double a10 = -100/ (double) 13*a00 + 5273/ (double) 1040;
	const double a01 = 3;
	const double a20 = 4;
	const double a11 = 2;
	const double a02 = 4;
	const double default_x = 0.13;
	const double default_y = 0.15;

	{
		const double x = default_x;
		const double disc_disc = 16*std::pow(a02, 2)*std::pow(a10, 2) - 16*a01*a02*a10*a11 + 16*a00*a02*std::pow(a11, 2) + 16*std::pow(a01, 2)*a02*a20 - 64*a00*std::pow(a02, 2)*a20;
		const double disc = (std::pow(a11, 2) - 4*a02*a20)*std::pow(x, 2) + std::pow(a01, 2) - 4*a00*a02 - 2*(2*a02*a10 - a01*a11)*x;
		REQUIRE(std::abs(disc_disc) < LARGER_TOL);
		REQUIRE(std::abs(disc + 1) < LARGER_TOL);
	}

	billiards::test::SolutionsTracker2d tracker;
	solve_22(a00, a10, a01, a20, a11, a02, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 1);
	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(
			a00 +
			a10 * p.first + a01 * p.second +
			a20 * p.first * p.first + a11 * p.first * p.second + a02 * p.second * p.second) < LARGER_TOL);
	}
}


TEST_CASE("system 22 discriminant has two roots and default_x does not work", "[solve_22]") {
	/*
ass = {'a01': 3, 'a20': 4, 'a11': 2, 'a02': 4, 'x': 0.13}
a10_sol = solve(disc.substitute(**ass) == -1, a10)[0].right()
a00_sol = solve(disc_disc.substitute(**ass, a10=a10_sol) == 1, a00)[0].right()
print('\tconst double a00 = ' + convert_powers(a00_sol) + ';')
print('\tconst double a10 = ' + convert_powers(a10_sol) + ';')

*/
	const double a00 = -13/ (double) 1600*sqrt(241) + 5507/ (double) 8000;
	const double a10 = -100/ (double) 13*a00 + 5273/ (double) 1040;
	const double a01 = 3;
	const double a20 = 4;
	const double a11 = 2;
	const double a02 = 4;
	const double default_x = 0.13;
	const double default_y = 0.15;

	{
		const double x = default_x;
		const double disc_disc = 16*std::pow(a02, 2)*std::pow(a10, 2) - 16*a01*a02*a10*a11 + 16*a00*a02*std::pow(a11, 2) + 16*std::pow(a01, 2)*a02*a20 - 64*a00*std::pow(a02, 2)*a20;
		const double disc = (std::pow(a11, 2) - 4*a02*a20)*std::pow(x, 2) + std::pow(a01, 2) - 4*a00*a02 - 2*(2*a02*a10 - a01*a11)*x;
		REQUIRE(std::abs(disc_disc - 1) < LARGER_TOL);
		REQUIRE(std::abs(disc + 1) < LARGER_TOL);
	}

	billiards::test::SolutionsTracker2d tracker;
	solve_22(a00, a10, a01, a20, a11, a02, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 2);
	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(
			a00 +
			a10 * p.first + a01 * p.second +
			a20 * p.first * p.first + a11 * p.first * p.second + a02 * p.second * p.second) < LARGER_TOL);
	}
}




///////////////////////////////////////////////////////
// solve_22
///////////////////////////////////////////////////////
TEST_CASE("system 22_22 happy case", "[solve_22_22]") {
	const double a00 = 1;
	const double a10 = -2;
	const double a01 = 3;
	const double a20 = -4;
	const double a11 = 5;
	const double a02 = -6;

	const double b00 = 6;
	const double b10 = -5;
	const double b01 = 4;
	const double b20 = -3;
	const double b11 = 2;
	const double b02 = -1;

	const double default_x = 0.13;
	const double default_y = 0.15;

	billiards::test::SolutionsTracker2d tracker;
	solve_22_22(
		a00, a10, a01, a20, a11, a02,
		b00, b10, b01, b20, b11, b02,
		default_x, default_y, tracker);
//	REQUIRE(tracker.count() == 4);
	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(
			a00 +
			a10 * p.first + a01 * p.second +
			a20 * p.first * p.first + a11 * p.first * p.second + a02 * p.second * p.second) < LARGER_TOL);
		REQUIRE(std::abs(
			b00 +
			b10 * p.first + b01 * p.second +
			b20 * p.first * p.first + b11 * p.first * p.second + b02 * p.second * p.second) < LARGER_TOL);
	}
}


void intersect_circles(
	const double c1x, const double c1y, const double r1, const double ax1, const double ax2,
	const double default_x, const double default_y,
	const billiards::test::SolutionsTracker2d& tracker
) {
	/*
x, y = var('x y')
cx, cy, r0, ax1, ax2 = var('cx cy r ax1 ax2')
((x - cx)^2 / ax1^2 + (y - cy)^2 / ax2^2 - r^2).expand()
# cx^2 + cy^2 - r^2 - 2*cx*x + x^2 - 2*cy*y + y^2
	 */
	const double c0x = 0;
	const double c0y = 0;
	const double r0 = 1;

	const double b00 = c0x * c0x + c0y * c0y - r0 * r0;
	const double b10 = -2 * c0x;
	const double b01 = -2 * c0y;
	const double b20 = 1;
	const double b11 = 0;
	const double b02 = 1;

	const double ax12 = ax1 * ax1;
	const double ax22 = ax2 * ax2;
	const double a00 = c1x * c1x / ax12 + c1y * c1y / ax22 - r1 * r1;
	const double a10 = -2 * c1x / ax12;
	const double a01 = -2 * c1y / ax22;
	const double a20 = 1 / ax12;
	const double a11 = 0;
	const double a02 = 1 / ax22;

//	const double default_x = 0.13;
//	const double default_y = 0.15;

	solve_22_22(
		a00, a10, a01, a20, a11, a02,
		b00, b10, b01, b20, b11, b02,
		default_x, default_y, tracker);

	for (const std::pair<double, double>& p : *tracker.sols) {
		REQUIRE(std::abs(
			a00 +
			a10 * p.first + a01 * p.second +
			a20 * p.first * p.first + a11 * p.first * p.second + a02 * p.second * p.second) < LARGER_TOL);
		REQUIRE(std::abs(
			b00 +
			b10 * p.first + b01 * p.second +
			b20 * p.first * p.first + b11 * p.first * p.second + b02 * p.second * p.second) < LARGER_TOL);
	}
}

TEST_CASE("system 22_22 with circles no solutions", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;
	billiards::test::SolutionsTracker2d tracker;
	intersect_circles(3, 0, 1, 1, 1, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 0);
}

TEST_CASE("system 22_22 with circles one solution", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;
	billiards::test::SolutionsTracker2d tracker;
	intersect_circles(2, 0, 1, 1, 1, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(std::pair<double, double>{1, 0}));
}

TEST_CASE("system 22_22 with circles two solutions", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;
	billiards::test::SolutionsTracker2d tracker;
	intersect_circles(1, 0, 1, 1, 1, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 2);
	REQUIRE(tracker.contains(std::pair<double, double>{0.5, +std::sqrt(3) / 2.0}));
	REQUIRE(tracker.contains(std::pair<double, double>{0.5, -std::sqrt(3) / 2.0}));
}

TEST_CASE("system 22_22 with circles four solutions", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;
	billiards::test::SolutionsTracker2d tracker;
	intersect_circles(0, 0, 1, 2, 0.5, default_x, default_y, tracker);
	REQUIRE(tracker.count() == 4);
//	REQUIRE(tracker.contains(std::pair<double, double>{0.5, +std::sqrt(3) / 2.0}));
//	REQUIRE(tracker.contains(std::pair<double, double>{0.5, -std::sqrt(3) / 2.0}));
}

TEST_CASE("system 22_22 with circles with infinitely many solutions", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;
	billiards::test::SolutionsTracker2d tracker;
	intersect_circles(0, 0, 1, 1, 1, default_x, default_y, tracker);
	REQUIRE(tracker.count() >= 1);
	REQUIRE((
		tracker.contains(std::pair<double, double>{default_x, std::sqrt(1 - default_x * default_x)}) ||
		tracker.contains(std::pair<double, double>{std::sqrt(1 - default_y * default_y), default_y})));
}

TEST_CASE("system 22_22 with circles with infinitely many solutions, no default value", "[solve_22_22]") {
	const double default_x = 30;
	const double default_y = 20;
	billiards::test::SolutionsTracker2d tracker;
	intersect_circles(0, 0, 1, 1, 1, default_x, default_y, tracker);
	REQUIRE(tracker.count() >= 1);
}


TEST_CASE("system 22_22 all zeros", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;

	billiards::test::SolutionsTracker2d tracker;
	solve_22_22(
		0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		default_x, default_y, tracker);
	REQUIRE(tracker.count() >= 1);
	REQUIRE(tracker.contains(std::pair<double, double>{default_x, default_y}));
}

TEST_CASE("system 22_22 almost all zeros but a constant term", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;

	billiards::test::SolutionsTracker2d tracker;
	solve_22_22(
		1, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0,
		default_x, default_y, tracker);
	REQUIRE(tracker.count() == 0);
}

TEST_CASE("system 22_22 almost all zeros, but a quadratic term", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;

	billiards::test::SolutionsTracker2d tracker;
	solve_22_22(
		0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 0, 1,
		default_x, default_y, tracker);
	REQUIRE(tracker.count() >= 1);
}

TEST_CASE("system 22_22 almost all zeros, but a quadratic cross term", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;

	billiards::test::SolutionsTracker2d tracker;
	solve_22_22(
		0, 0, 0, 0, 0, 0,
		1, 0, 0, 0, 1, 0,
		default_x, default_y, tracker);
	REQUIRE(tracker.count() >= 1);
}
TEST_CASE("system 22_22 linear and quadratic", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;

	billiards::test::SolutionsTracker2d tracker;
	solve_22_22(
		-25, 0, 0, 1, 0, 1,
		-3, 1, 0, 0, 0, 0,
		default_x, default_y, tracker);
	REQUIRE(tracker.count() >= 2);
	REQUIRE(tracker.contains(std::pair<double, double>{3, 4}));
	REQUIRE(tracker.contains(std::pair<double, double>{3, -4}));
}
TEST_CASE("system 22_22 hyperbola and sphere", "[solve_22_22]") {
	const double default_x = 0.13;
	const double default_y = 0.15;

	billiards::test::SolutionsTracker2d tracker;
	solve_22_22(
		-4, 0, 0, 1, 0, 1,
		-1, 0, 0, -1, 0, 1,
		default_x, default_y, tracker);
	REQUIRE(tracker.count() >= 4);
}