cmake_minimum_required(VERSION 3.20)
project(idea)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE on)


add_subdirectory(billiards-common)
add_subdirectory(billiards-layouts-api)
add_subdirectory(billiards-projection-api)
add_subdirectory(billiards-shots-api)

# TODO: rename shots and gphx
# TODO: use Dimensions....
# TODO: make color use uint8_t