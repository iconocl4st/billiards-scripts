//
// Created by thallock on 10/3/21.
//

#include "catch2/catch.hpp"

// Define this in only one file to add 'main'
//#define CATCH_CONFIG_MAIN

#include "billiards_common/math/polynomial.h"
#include "SolutionsTracker.h"

using namespace billiards::math;


////////////////////////////////////////////////////////
// Linear
////////////////////////////////////////////////////////
TEST_CASE("Solve linear happy", "[solve_1]") {
	const billiards::test::SolutionsTracker tracker;
	solve_1n(3.0, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(std::abs(3.0 + tracker[0]) < LARGER_TOL);
}
TEST_CASE("Solve linear default", "[solve_1]") {
	const billiards::test::SolutionsTracker tracker;
	solve_1(0.0, 0.0, 13.0, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(std::abs(-13 + tracker[0]) < LARGER_TOL);
}
TEST_CASE("Solve linear no solutions", "[solve_1]") {
	const billiards::test::SolutionsTracker tracker;
	solve_1(1.0, 0.0, 13.0, tracker);
	REQUIRE(tracker.count() == 0);
}


////////////////////////////////////////////////////////
// Quadratic
////////////////////////////////////////////////////////
TEST_CASE("Solve quadratic happy", "[solve_2]") {
	const double c = 3;
	const double r1 = 3;
	const double r2 = -2;
	const double c0 = c * r1 * r2;
	const double c1 = -c * (r1 + r2);
	const double c2 = c;
	const billiards::test::SolutionsTracker tracker;
	solve_2(c0, c1, c2, 13, tracker);
	REQUIRE(tracker.count() == 2);
	REQUIRE(tracker.contains(r1));
	REQUIRE(tracker.contains(r2));
}
TEST_CASE("Solve quadratic default", "[solve_2]") {
	const billiards::test::SolutionsTracker tracker;
	solve_2(0.0, 0.0, 0.0, 13.0, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(13.0));
}
TEST_CASE("Solve quadratic no solutions", "[solve_2]") {
	const billiards::test::SolutionsTracker tracker;
	solve_2(1.0, 0.0, 1.0, 13.0, tracker);
	REQUIRE(tracker.count() == 0);
}
TEST_CASE("Solve linear quadratic", "[solve_2]") {
	const billiards::test::SolutionsTracker tracker;
	solve_2(1.0, 3.0, 0.0, 13.0, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE((1.0 + 3.0 * tracker[0]) < LARGER_TOL);
}
TEST_CASE("Solve quadratic with 0", "[solve_2]") {
	const billiards::test::SolutionsTracker tracker;
	solve_2(0.0, 3.0, 1.0, 13.0, tracker);
	REQUIRE(tracker.count() == 2);
	REQUIRE((0.0 + 3.0 * tracker[0] + 1.0 * tracker[0] * tracker[0]) < LARGER_TOL);
	REQUIRE((0.0 + 3.0 * tracker[1] + 1.0 * tracker[0] * tracker[0]) < LARGER_TOL);
}
TEST_CASE("Solve quadratic repeated solutions", "[solve_2]") {
	// TODO: Fails at 0
	const double c = 3;
	const double r1 = 3;
	const double r2 = 3;
	const double c0 = c * r1 * r2;
	const double c1 = -c * (r1 + r2);
	const double c2 = c;
	const billiards::test::SolutionsTracker tracker;
	solve_2(c0, c1, c2, 13, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(r1));
}


////////////////////////////////////////////////////////
// Cubic
////////////////////////////////////////////////////////
TEST_CASE("Solve cubic happy", "[solve_3]") {
	const double c = 3;
	const double r1 = 4;
	const double r2 = -2;
	const double r3 = 5;
	const double c0 = -c * r1 * r2 * r3;
	const double c1 = c * (r1 * r2 + r1 * r3 + r2 * r3);
	const double c2 = -c * (r1 + r2 + r3);
	const double c3 = c;
	const billiards::test::SolutionsTracker tracker;
	solve_3(c0, c1, c2, c3, 13, tracker);
	REQUIRE(tracker.count() == 3);
	REQUIRE(tracker.contains(r1));
	REQUIRE(tracker.contains(r2));
	REQUIRE(tracker.contains(r3));
}
TEST_CASE("Solve cubic default", "[solve_3]") {
	const billiards::test::SolutionsTracker tracker;
	solve_3(0.0, 0.0, 0.0, 0.0, 13.0, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(13.0));
}
TEST_CASE("Solve cubic one solutions", "[solve_3]") {
	const billiards::test::SolutionsTracker tracker;
	solve_3(-5, 0, 0, 1, 13, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(std::cbrt(-5)));
}
TEST_CASE("Solve quadratic cubic", "[solve_3]") {
	const billiards::test::SolutionsTracker tracker;
	solve_3(-4, 0, 1, 0, 13, tracker);
	REQUIRE(tracker.count() == 2);
	REQUIRE(tracker.contains(2));
	REQUIRE(tracker.contains(-2));
}
TEST_CASE("Solve linear cubic", "[solve_3]") {
	// TODO: fails at 0
	const billiards::test::SolutionsTracker tracker;
	solve_3(-4, 1, 0, 0, 13, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(4));
}
TEST_CASE("Solve cubic repeated solutions", "[solve_3]") {
	// TODO: Fails at 0
	const double c = 3;
	const double r1 = 3;
	const double r2 = 3;
	const double r3 = 5;
	const double c0 = -c * r1 * r2 * r3;
	const double c1 = c * (r1 * r2 + r1 * r3 + r2 * r3);
	const double c2 = -c * (r1 + r2 + r3);
	const double c3 = c;
	const billiards::test::SolutionsTracker tracker;
	solve_3(c0, c1, c2, c3, 13, tracker);
	// TODO: Should be 2:
	REQUIRE(tracker.count() == 3);
	REQUIRE(tracker.contains(r1));
	REQUIRE(tracker.contains(r3));
}

TEST_CASE("Solve cubic all repeated solutions", "[solve_3]") {
	// TODO: Fails at 0
	const double c = 3;
	const double r1 = 3;
	const double r2 = 3;
	const double r3 = 3;
	const double c0 = -c * r1 * r2 * r3;
	const double c1 = c * (r1 * r2 + r1 * r3 + r2 * r3);
	const double c2 = -c * (r1 + r2 + r3);
	const double c3 = c;
	const billiards::test::SolutionsTracker tracker;
	solve_3(c0, c1, c2, c3, 13, tracker);
	// TODO: Should be 2:
	REQUIRE(tracker.count() == 3);
	REQUIRE(tracker.contains(r1));
	for (int i = 0; i < 3; i++) {
		REQUIRE(std::abs(c0 + c1 * tracker[i] + c2 * tracker[i] * tracker[i] + c3 * tracker[i] * tracker[i] * tracker[i]) < LARGER_TOL);
	}
}


//////////////////////////////////////////////////////////
//// Quartic
//////////////////////////////////////////////////////////
TEST_CASE("Solve quartic happy", "[solve_4]") {
	const double c = 3;
	const double r1 = 4;
	const double r2 = -2;
	const double r3 = 5;
	const double r4 = -17;
	const double c0 = c * r1 * r2 * r3 * r4;
	const double c1 = -c * (r1 * r2 * r3 + r1 * r2 * r4 + r1 * r3 * r4 + r2 * r3 * r4);
	const double c2 = c * (r1 * r2 + r1 * r3 + r1 * r4 + r2 * r3 + r2 * r4 + r3 * r4);
	const double c3 = -c * (r1 + r2 + r3 + r4);
	const double c4 = c;
	const billiards::test::SolutionsTracker tracker;
	solve_4(c0, c1, c2, c3, c4, 13, tracker);
	REQUIRE(tracker.count() == 4);
	REQUIRE(tracker.contains(r1));
	REQUIRE(tracker.contains(r2));
	REQUIRE(tracker.contains(r3));
	REQUIRE(tracker.contains(r4));
}
TEST_CASE("Solve quartic default", "[solve_4]") {
	const billiards::test::SolutionsTracker tracker;
	solve_4(0.0, 0.0, 0.0, 0.0, 0.0, 13.0, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(13.0));
}
TEST_CASE("Solve quartic no solutions", "[solve_4]") {
	const billiards::test::SolutionsTracker tracker;
	solve_4(1, 0, 0, 0, 1, 13, tracker);
	REQUIRE(tracker.count() == 0);
}
TEST_CASE("Solve linear quartic", "[solve_4]") {
	const billiards::test::SolutionsTracker tracker;
	solve_4(1, 1, 0, 0, 0, 13, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(-1));
}
TEST_CASE("Solve quadratic quartic", "[solve_4]") {
	const billiards::test::SolutionsTracker tracker;
	solve_4(-16, 0, 1, 0, 0, 13, tracker);
	REQUIRE(tracker.count() == 2);
	REQUIRE(tracker.contains(4));
	REQUIRE(tracker.contains(-4));
}
TEST_CASE("Solve cubic quartic", "[solve_4]") {
	const billiards::test::SolutionsTracker tracker;
	solve_4(27, 0, 0, 1, 0, 13, tracker);
	REQUIRE(tracker.count() == 1);
	REQUIRE(tracker.contains(3));
}
TEST_CASE("Solve root repeated twice", "[solve_4]") {
	// TODO: Fails at 0
	const double c = 3;
	const double r1 = 4;
	const double r2 = 4;
	const double r3 = 5;
	const double r4 = 5;
	const double c0 = c * r1 * r2 * r3 * r4;
	const double c1 = -c * (r1 * r2 * r3 + r1 * r2 * r4 + r1 * r3 * r4 + r2 * r3 * r4);
	const double c2 = c * (r1 * r2 + r1 * r3 + r1 * r4 + r2 * r3 + r2 * r4 + r3 * r4);
	const double c3 = -c * (r1 + r2 + r3 + r4);
	const double c4 = c;
	const billiards::test::SolutionsTracker tracker;
	solve_4(c0, c1, c2, c3, c4, 13, tracker);
	for (const double s : *tracker.sols) {
		REQUIRE(std::abs(c0 + c1 * s + c2 * s * s + c3 * s * s * s + c4 * s * s * s * s) < LARGER_TOL);
	}
	REQUIRE(tracker.count() >= 2);
	REQUIRE(tracker.contains(r1));
	REQUIRE(tracker.contains(r4));
}
TEST_CASE("Solve root repeated 3 times", "[solve_4]") {
	// TODO: Fails at 0
	const double c = 3;
	const double r1 = 4;
	const double r2 = 4;
	const double r3 = 4;
	const double r4 = 5;
	const double c0 = c * r1 * r2 * r3 * r4;
	const double c1 = -c * (r1 * r2 * r3 + r1 * r2 * r4 + r1 * r3 * r4 + r2 * r3 * r4);
	const double c2 = c * (r1 * r2 + r1 * r3 + r1 * r4 + r2 * r3 + r2 * r4 + r3 * r4);
	const double c3 = -c * (r1 + r2 + r3 + r4);
	const double c4 = c;
	const billiards::test::SolutionsTracker tracker;
	solve_4(c0, c1, c2, c3, c4, 13, tracker);
	for (const double s : *tracker.sols) {
		REQUIRE(std::abs(c0 + c1 * s + c2 * s * s + c3 * s * s * s + c4 * s * s * s * s) < LARGER_TOL);
	}
	REQUIRE(tracker.count() >= 2);
	REQUIRE(tracker.contains(r1));
	REQUIRE(tracker.contains(r4));
}
TEST_CASE("Solve all repeated solutions", "[solve_4]") {
	// TODO: Fails at 0
	const double c = 3;
	const double r1 = 4;
	const double r2 = 4;
	const double r3 = 4;
	const double r4 = 4;
	const double c0 = c * r1*r2*r3*r4;
	const double c1 = -c * (r1*r2*r3 + r1*r2*r4 + r1*r3*r4 + r2*r3*r4);
	const double c2 = c * (r1*r2 + r1*r3 + r2*r3 + r1*r4 + r2*r4 + r3*r4);
	const double c3 = -c * (r1 + r2 + r3 + r4);
	const double c4 = c;
	const double val = c0 + c1 * 4  + c2 * 4 * 4 + c3 * 4 * 4 * 4 + c4 * 4 * 4 * 4 * 4;
	std::cout << val << std::endl;
	const billiards::test::SolutionsTracker tracker;
	solve_4(c0, c1, c2, c3, c4, 13, tracker);
	for (const double s : *tracker.sols) {
		REQUIRE(std::abs(c0 + c1 * s + c2 * s * s + c3 * s * s * s + c4 * s * s * s * s) < LARGER_TOL);
	}
	REQUIRE(tracker.count() >= 1);
	REQUIRE(tracker.contains(r1));
}

TEST_CASE("Cubes", "[complex]") {
	billiards::math::Complex c{3.4, -2.9};
	billiards::math::Complex r = c.pow(1 / 3.0);
	REQUIRE(std::abs(((r * r * r) - c).r()) < LARGER_TOL);
}

//TEST_CASE("Solve quartic with two solutions", "[solve_4]") {
//TEST_CASE("Solve quartic with three solutions", "[solve_4]") {

//x = var('x')
//expand((x - 3) * (x + 2) * (x-7) * (x + 4))
//expand((x - 3) * (x + 2) * (x-7))
//
// */
//		const double c4 = 1;
//		const double c3 = -4;
//		const double c2 = -31;
//		const double c1 = 46;
//		const double c0 = 168;
//
//		{
//			for (const double x: std::array<double, 4>{3, -2, 7, -4}) {
//				std::cout << (c0 + c1 * x + c2 * x * x + c3 * x * x * x + c4 * x * x * x * x) << std::endl;
//			}
//			std::cout << "-" << std::endl;
//			for (const double x: std::array<double, 3>{3, -2, 7}) {
//				std::cout << (42 + x - 8 * x * x + x * x * x) << std::endl;
//			}
//			std::cout << "-" << std::endl;
//		}
//
//		const double P = -6;
//		const double Q = 4;
//		solve_3_simp(P, Q, [P, Q](double x) {
//			std::cout << "Simplified cubic solution: " << x << std::endl;
//			std::cout << "value: " << (x * x * x + P * x - Q) << std::endl;
//		});
//
//		solve_2(
//			12, -7, 1.0,
//			0,
//			[](double x) {
//				std::cout << "Quadratic solution: " << x << std::endl;
//				std::cout << "p(x) = " << (12 - 7 * x + x * x) << std::endl;
//			});
//
//		solve_3(
//			42, 1.0, -8, 1.0,
//			0,
//			[](double x) {
//				std::cout << "Cubic solution: " << x << std::endl;
//				std::cout << "p(x) = " << (42 + x - 8 * x * x + x * x * x) << std::endl;
//			});
//
//		solve_4(
//			c0, c1, c2, c3, c4,
//			0,
//			[c0, c1, c2, c3, c4](double x) {
//				std::cout << "Quartic solution: " << x << std::endl;
//				std::cout << "p(x) = " << (c0 + c1 * x + c2 * x * x + c3 * x * x * x + c4 * x * x * x * x) << std::endl;
//			});
//	}
//
//	if (0)
//	{
//		const double c0 = -0.034133333333333335;
//		const double c1 = 0.17066666666666666;
//		const double c2 = 0;
//		const double c3 = -0.53333333333333333;
//		// -0.034133333333333335 + 0.17066666666666666 * x + -0.53333333333333333 * x^2 + x^4
//
//		solve_4n(
//			c0, c1, c2, c3,
//			[c0, c1, c2, c3](double x) {
//				std::cout << "Quartic solution: " << x << std::endl;
//				std::cout << "p(x) = " << (c0 + c1 * x + c2 * x * x + c3 * x * x * x + x * x * x * x) << std::endl;
//			});
//	}
//
//	{
//		const double a00 = -1;
//		const double a10 = 0;
//		const double a01 = 0;
//		const double a20 = 1;
//		const double a11 = 0;
//		const double a02 = 1;
//
//		const double b00 = -1;
//		const double b10 = 0;
//		const double b01 = 0;
//		const double b20 = 1;
//		const double b11 = 0;
//		const double b02 = 0;
//		solve_22(
//			a00, a10, a01, a20, a11, a02,
//			b00, b10, b01, b20, b11, b02,
//			std::nan(""), std::nan(""),
//			[
//				a00, a10, a01, a20, a11, a02,
//				b00, b10, b01, b20, b11, b02](const double x, const double y
//			) {
//				std::cout << "=========================================================" << std::endl;
//				std::cout << "systems solution:" << std::endl;
//				std::cout << x << ", " << y << std::endl;
//
//				std::cout << "eq1: " << a00 + a10 * x + a01 * y + a20 * x * x + a11 * x * y + a02 * y * y << std::endl;
//				std::cout << "eq2: " << b00 + b10 * x + b01 * y + b20 * x * x + b11 * x * y + b02 * y * y << std::endl;
//				std::cout << "=========================================================" << std::endl;
//			}
//		);
//
//		{
//			const double x = 1;
//			const double y = 0;
//			std::cout << "=========================================================" << std::endl;
//			std::cout << "some solution" << std::endl;
//			std::cout << "eq1: " << a00 + a10 * x + a01 * y + a20 * x * x + a11 * x * y + a02 * y * y << std::endl;
//			std::cout << "eq2: " << b00 + b10 * x + b01 * y + b20 * x * x + b11 * x * y + b02 * y * y << std::endl;
//			std::cout << "=========================================================" << std::endl;
//		}
//	}
//
//	return 0;
//}
