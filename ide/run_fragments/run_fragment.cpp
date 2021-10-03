//
// Created by thallock on 10/1/21.
//

#include <iostream>

#include "./run_glance_code.h"
#include "./run_orient_pocket.h"
#include "./run_bank_shot.h"

int main(int argc, char **argv) {
	std::cout << "Running a fragment" << std::endl;

//	billiards::fragments::run_bank_code();
//	billiards::fragments::run_glance_code();
	billiards::fragments::run_orient_pocket();
	return 0;
}
