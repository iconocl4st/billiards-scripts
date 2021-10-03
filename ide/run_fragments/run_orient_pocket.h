//
// Created by thallock on 10/2/21.
//

#ifndef IDEA_RUN_ORIENT_POCKET_H
#define IDEA_RUN_ORIENT_POCKET_H

#include "math/orient_pocket.h"
#include "billiards_common/geometry/geometry.h"

namespace billiards::fragments {


	void run_orient_pocket() {

		const geometry::MaybePoint source{20, 20};
		const geometry::MaybePoint segment1{0, 0};
		const geometry::MaybePoint segment2{50, 0};
		const geometry::MaybeDouble radius{10};
		const geometry::MaybePoint solution = shots::math::orient_pocket(
			source, segment1, segment2, radius);

		if (!solution.is_valid()) {
			std::cout << "No solution" << std::endl;
			return;
		}

		std::cout << "Solution: " << solution << std::endl;
		const auto travel_line = geometry::through(solution, source);
		const auto normal_line = geometry::orthogonal_at(travel_line, segment1);
		const auto intersection = geometry::intersection(normal_line, travel_line);
		std::cout << "Distance: " << (intersection - segment1).norm() << std::endl;
	}
}
#endif //IDEA_RUN_ORIENT_POCKET_H
