//
// Created by thallock on 10/29/21.
//


#include "catch2/catch.hpp"

#include "billiards_common/utils/fs.h"

using namespace billiards::utils;

TEST_CASE("ensure directory of file", "[ensure_directory]") {
    const auto b = ensure_directory("/work/pool/repos/billiards-scripts/Makefile");
    REQUIRE(!b);
}

TEST_CASE("ensure directory", "[ensure_directory]") {
    const auto b = ensure_directory("/work/pool/repos/billiards-scripts/tmp_dir");
    REQUIRE(b);
}

TEST_CASE("list directory", "[list_directory]") {
    int count = 0;
    const auto rec = [&](const std::string& p) {
        std::cout << p << std::endl;
        count++;
    };
    list_directory("/work/pool/repos/billiards-scripts", rec);

    REQUIRE(count > 0);
}
