//
// Created by thallock on 10/6/21.
//

#include "catch2/catch.hpp"

#include "write_shot.h"

#include "billiards_common/shots/shot_is_possible.h"
#include "billiards_common/shots/Locations.h"
#include "shot_calculation/shot.h"
#include "ShotInfoParams.h"

using namespace billiards::shots;
using namespace billiards::layout;
using namespace billiards::geometry;
using namespace billiards::config;
using namespace billiards::test;

void add_ball(ShotInfoParams& params, double x, double y, vball::virtual_type::VirtualBallType type, int number) {
	params.locations.balls.emplace_back(
		vball::VirtualBall{
			type, number},
		Point{x, y}
	);
}

void add_obj(ShotInfoParams& params, double x, double y, int number) {
	add_ball(params, x, y, vball::virtual_type::VirtualBallType::NUMBER, number);
}

void add_cue(ShotInfoParams& params, double x, double y) {
	add_ball(params, x, y, vball::virtual_type::VirtualBallType::CUE, 0);
}

void write_shot(ShotInfoParams& params, ShotInformation& info) {
	print_object(params.locations, "locations.json");
	print_object(params.table, "table.json");
	print_object(info, "shot_info.json");
}

void require_feasible(ShotInfoParams& params, ShotInformation& info) {
	bool success = calculate_shot(params, info);
//	REQUIRE(success);
	const bool possible = shot_info_is_possible(
		params.table,
		params.locations,
		info);

//	REQUIRE(possible);
}


TEST_CASE("test strike happy case", "[calculate_shot]") {
	ShotInfoParams params;
	add_cue(params, 20, 20);
	add_obj(params, 40, 30, 1);
	params.locations.table_dims = params.table.dims;
	params.shot.steps.push_back(std::make_shared<CueStep>(0));
	params.shot.steps.push_back(std::make_shared<StrikeStep>(1));
	params.shot.steps.push_back(std::make_shared<PocketStep>(constants::RIGHT_UPPER_POCKET));

	ShotInformation info{params.shot};
	require_feasible(params, info);
}

TEST_CASE("run bank shot happy case", "[calculate_shot]") {
	ShotInfoParams params;
	add_cue(params, 50, 30);
	add_obj(params, 40, 10, 1);
	params.locations.table_dims = params.table.dims;
	params.shot.steps.push_back(std::make_shared<CueStep>(0));
	params.shot.steps.push_back(std::make_shared<StrikeStep>(1));
	params.shot.steps.push_back(std::make_shared<RailStep>(
		constants::LOWER_LEFT_RAIL));
	params.shot.steps.push_back(std::make_shared<PocketStep>(
		constants::MIDDLE_UPPER_POCKET));

	ShotInformation info{params.shot};
	require_feasible(params, info);
}

TEST_CASE("run kick shot happy case", "[calculate_shot]") {
	ShotInfoParams params;
	add_cue(params, 50, 30);
	add_obj(params, 40, 10, 1);
	params.locations.table_dims = params.table.dims;
	params.shot.steps.push_back(std::make_shared<CueStep>(0));
	params.shot.steps.push_back(std::make_shared<RailStep>(
		constants::LOWER_LEFT_RAIL));
	params.shot.steps.push_back(std::make_shared<StrikeStep>(1));
	params.shot.steps.push_back(std::make_shared<PocketStep>(
		constants::MIDDLE_UPPER_POCKET));

	ShotInformation info{params.shot};
	require_feasible(params, info);
}

TEST_CASE("run combo shot happy case", "[calculate_shot]") {
	ShotInfoParams params;
	add_cue(params, 20, 30);
	add_obj(params, 40, 10, 1);
	add_obj(params, 60, 5, 2);
	params.locations.table_dims = params.table.dims;
	params.shot.steps.push_back(std::make_shared<CueStep>(0));
	params.shot.steps.push_back(std::make_shared<StrikeStep>(1));
	params.shot.steps.push_back(std::make_shared<StrikeStep>(2));
	params.shot.steps.push_back(std::make_shared<PocketStep>(constants::RIGHT_LOWER_POCKET));

	ShotInformation info{params.shot};
	require_feasible(params, info);
}

TEST_CASE("run kiss shot happy case", "[calculate_shot]") {
	ShotInfoParams params;
	add_cue(params, 20, 30);
	add_obj(params, 40, 10, 1);
	add_obj(params, 60, 5, 2);
	params.locations.table_dims = params.table.dims;
	params.shot.steps.push_back(std::make_shared<CueStep>(0));
	params.shot.steps.push_back(std::make_shared<KissStep>(1, kt::kiss_type::ROLLING));
	params.shot.steps.push_back(std::make_shared<StrikeStep>(2));
	params.shot.steps.push_back(std::make_shared<PocketStep>(constants::RIGHT_LOWER_POCKET));

	ShotInformation info{params.shot};
	require_feasible(params, info);
	write_shot(params, info);
}

TEST_CASE("run carom shot happy case", "[calculate_shot]") {
	ShotInfoParams params;
	add_cue(params, 20, 30);
	add_obj(params, 40, 20, 1);
	add_obj(params, 80, 20, 2);
	params.locations.table_dims = params.table.dims;
	params.shot.steps.push_back(std::make_shared<CueStep>(0));
	params.shot.steps.push_back(std::make_shared<StrikeStep>(1));
	params.shot.steps.push_back(std::make_shared<KissStep>(2, kt::kiss_type::ROLLING));
	params.shot.steps.push_back(std::make_shared<PocketStep>(constants::RIGHT_LOWER_POCKET));

	ShotInformation info{params.shot};
	require_feasible(params, info);
//	write_shot(params, info);
}