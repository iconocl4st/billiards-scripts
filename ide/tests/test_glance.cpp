//
// Created by thallock on 10/5/21.
//

#include "catch2/catch.hpp"

// Define this in only one file to add 'main'
//#define CATCH_CONFIG_MAIN

#include "math/rolling_glance.h"
#include "SolutionsTracker.h"

#define LARGER_TOLERANCE 1e-5


TEST_CASE("Solve get_glance_location happy case", "[get_glance_location]") {
	const double ax = -1;
	const double ay = 1;
	const double dx = 2;
	const double dy = -1;
	const double r = 0.5;

	billiards::shots::math::get_glance_location(
		ax, ay,
		dx, dy,
		r,
		[&](const billiards::shots::RollingGlanceCalculation& c) {
			const double x = c.loc.x;
			const double y = c.loc.y;

			std::cout << "x=" << x << ", y=" << y << std::endl;

			const double tx = y;
			const double ty = -x;
			const double s = (ax * tx + ay * ty) / (ax * ax + ay * ay);
			const double alpha = 2/7.0;

			REQUIRE(std::abs(x * x + y * y - r * r) < LARGER_TOLERANCE);
			REQUIRE(std::abs((s * ax) * (s * ax - tx) + (s * ay) * (s * ay - ty)) < LARGER_TOLERANCE);

			// These are checked within the method:
//			eq_dx = t * dx == (1 - alpha) * tx + alpha * s * ax
//			eq_dy = t * dy == (1 - alpha) * ty + alpha * s * ay
		}
	);
//	REQUIRE(tracker.count() == 1);
//	REQUIRE(std::abs(3.0 + tracker[0]) < LARGER_TOL);
}

