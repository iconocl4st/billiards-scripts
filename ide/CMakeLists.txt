cmake_minimum_required(VERSION 3.18)
project(idea)

include(FindPkgConfig)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE on)

add_subdirectory(billiards-common)
add_subdirectory(billiards-layouts-api)
add_subdirectory(billiards-projection-api)
add_subdirectory(billiards-shots-api)
add_subdirectory(billiards-graphics-api)

find_package(Boost COMPONENTS system REQUIRED)
find_package(Catch2 2 REQUIRED hints "/usr/lib64/cmake/Catch2/")
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

pkg_search_module(gmpxx REQUIRED gmpxx)


if (NOT DEFINED $ENV{REPOS})
	set($ENV{REPOS} "/work/pool/repos")
endif()

add_executable(
        run_tests

		${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tests/test_sparse_reduction.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tests/test_glance_system.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tests/test_gmp_poly.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tests/test_parse_poly.cpp
		${CMAKE_CURRENT_SOURCE_DIR}/tests/test_higher_order_polys.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_polynomials.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_systems.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_glance.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_shot_calculation.cpp
        tests/test_mpfr.cpp)
target_link_libraries(run_tests PRIVATE Eigen3::Eigen gmpxx gmp mpfr)
target_include_directories(
        run_tests

        PUBLIC

        ${PROJECT_SOURCE_DIR}/src

        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-attempts-api/src
        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-common/src
        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-config-api/src
        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-graphics-api/src
        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-layouts/src
        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-projection-api/src
        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-shots-api/src

        $ENV{REPOS}/billiards-scripts/ide/billiards-attempts-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-common/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-config-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-graphics-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-layouts-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-projection-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-shots-api/src
)

add_executable(
        run_fragment

        ${CMAKE_CURRENT_SOURCE_DIR}/run_fragments/run_fragment.cpp
#        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-shots-api/src/math/orient_pocket_gen/gen_sol_0.c
#        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-shots-api/src/math/orient_pocket_gen/gen_sol_1.c
#        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-shots-api/src/math/orient_pocket_gen/gen_sol_2.c
#        ${CMAKE_CURRENT_SOURCE_DIR}/billiards-shots-api/src/math/orient_pocket_gen/gen_sol_3.c

)

target_include_directories(
        run_fragment

        PUBLIC

        ${PROJECT_SOURCE_DIR}/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-attempts-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-common/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-config-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-graphics-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-layouts/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-projection-api/src
        $ENV{REPOS}/billiards-scripts/ide/billiards-shots-api/src
)


# TODO: rename shots and gphx
# TODO: use Dimensions
# TODO: make color use uint8_t
# TODO: make consistent enums
# TODO: make consistent serialization from subclasses
# TODO: Maybe remove the destination?
# TODO: check if shot is possible
# TODO: Make necessary elements parsed by json
# TODO: Make the layout fit the table
# TODO: Place the pocket points in a nicer place
# TODO: Only have one ghost ball (the center)
# TODO: params vs table vs locations
# TODO: parallelize build
# TODO: Build could be faster if the crow server had a wrapper
# TODO: orient pocket has a bug where it widens the angle...
# TODO: tests
# TODO: Maybe go to eigen
