//
// Created by thallock on 12/22/21.
//


#include "catch2/catch.hpp"

#include "billiards_common/geometry/Point.h"
#include "construct_map.h"
#include "billiards_common/config/TableGeometry.h"

using namespace billiards::geometry;
using namespace billiards::config;
using namespace billiards::project;

TEST_CASE("Find line", "[smh]") {
	const Point v1{5, 10};
	const Point v2{3, -6};
	const Point o{0, 0};

	Hyperplane h;

	std::cout << h.orient(v1.x, v1.y, v2.x, v2.y, o.x, o.y) << std::endl;
	std::cout << h << std::endl;

	std::cout << h.at(v1.x, v1.y) << std::endl;
	std::cout << h.at(v2.x, v2.y) << std::endl;
	std::cout << h.at(o.x, o.y) << std::endl;
}


TEST_CASE("Construct map", "[construct_map]") {
	MappingInformation info;
	TriangulationMap map;
	bool success = construct_map(info, map);
	REQUIRE(success);
}

TEST_CASE("Apply mapping", "[construct_map]") {
	MappingInformation info;
	TriangulationMap map;
	bool success = construct_map(info, map);
	REQUIRE(success);
	Point in{3, 4};
	Point out = map.map(in);
	std::cout << out << std::endl;
}

TEST_CASE("Apply mapping to out of bounds point", "[construct_map]") {
	MappingInformation info;
	TriangulationMap map;
	bool success = construct_map(info, map);
	REQUIRE(success);
	Point in{-3, -4};
	Point out = map.map(in);
	std::cout << out << std::endl;
}