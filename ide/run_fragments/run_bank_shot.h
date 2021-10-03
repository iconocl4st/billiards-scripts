//
// Created by thallock on 10/2/21.
//

#ifndef IDEA_RUN_BANK_SHOT_H
#define IDEA_RUN_BANK_SHOT_H


#include "common_fragments.h"

#include "billiards_common/shots/ShotInformation.h"
#include "billiards_common/shots/shot_is_possible.h"

#include "ShotInfoParams.h"
#include "calculate_shot.h"

namespace billiards::fragments {

	inline
	void run_bank_code() {
		shots::ShotInfoParams params;
		params.locations.balls.emplace_back(
			layout::vball::VirtualBall{
				layout::vball::virtual_type::CUE, 0},
			geometry::Point{20.0, 20.0}
		);
		params.locations.balls.emplace_back(
			layout::vball::VirtualBall{
				layout::vball::virtual_type::NUMBER, 1},
			geometry::Point{40.0, 20.0}
		);
		params.locations.table_dims = params.table.dims;
		params.shot.steps.push_back(std::make_shared<shots::CueStep>(0));
		params.shot.steps.push_back(std::make_shared<shots::StrikeStep>(1));
		params.shot.steps.push_back(std::make_shared<shots::RailStep>(
			config::constants::LOWER_RIGHT_RAIL));
		params.shot.steps.push_back(std::make_shared<shots::PocketStep>(
			config::constants::RIGHT_UPPER_POCKET));

		shots::ShotInformation info{params.shot};
		shots::calculate_shot(params, info);

		print_object(params.locations, "locations.json");
		print_object(params.table, "table.json");
		print_object(info, "shot_info.json");

		const bool possible = shots::shot_info_is_possible(
			params.table,
			params.locations,
			info);
		std::cout << "Shot is possible: " << possible << std::endl;
	}
}



#endif //IDEA_RUN_BANK_SHOT_H
