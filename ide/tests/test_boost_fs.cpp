//
// Created by thallock on 10/29/21.
//


#include "catch2/catch.hpp"

#include <iostream>
#include <boost/filesystem.hpp>


bool ensure_directory(const std::string& path) {
    boost::filesystem::path dir{path};
    if (!boost::filesystem::exists(dir)) {
        boost::filesystem::create_directories(dir);
    }
    return boost::filesystem::is_directory(dir);
}

void list_directory(const std::string& path, const std::function<void(const std::string&)>& rec) {
    if (!ensure_directory(path)) {
        return;
    }
    boost::filesystem::path dir{path};
    for (auto& child : boost::filesystem::directory_iterator(dir)) {
        rec(child.path().string());
    }
}



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
