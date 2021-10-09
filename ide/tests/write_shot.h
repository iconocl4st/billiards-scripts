//
// Created by thallock on 10/6/21.
//

#ifndef IDEA_WRITE_SHOT_H
#define IDEA_WRITE_SHOT_H

#include <fstream>
#include <iostream>

#include "billiards_common/utils/dump.h"


namespace billiards::test {

	const std::string OUT_DIR = "/mnt/1f0ab4b3-c472-49e1-92d8-c0b5664f7fdb/ProjectsForFun/Pool/repos/billiards-scripts/scripts/node-app/inputs";

	inline
	void print_object(const json::Serializable& obj, const std::string& filename) {
		const std::string output_str = json::pretty_dump(obj);
		std::stringstream out_file_name;
		out_file_name << OUT_DIR << "/" << filename;
		std::ofstream out_file{out_file_name.str()};
		out_file << output_str;

		std::cout << "Writing to " << filename << ":" << std::endl;
//		std::cout << output_str << std::endl;
	}
}

#endif //IDEA_WRITE_SHOT_H
