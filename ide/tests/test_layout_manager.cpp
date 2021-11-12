//
// Created by thallock on 11/9/21.
//

#include "catch2/catch.hpp"

#include "billiards_common/unql/Store.h"

//using namespace billiards::layout;

class TestClass : public billiards::json::Serializable {
public:
	std::string str;

	explicit TestClass(std::string val) : str{std::move(val)} {}
	TestClass() : TestClass{"No value provided"} {}

	~TestClass() override = default;

	void parse(const nlohmann::json& value, billiards::json::ParseResult& status) override {
		ENSURE_STRING(status, value, "value", "Must have a string");
		str = value["value"].get<std::string>();
	}

	void to_json(billiards::json::SaxWriter& writer) const override {
		writer.begin_object();
		writer.field("value", str);
		writer.end_object();
	}
};

TEST_CASE("do some stuff", "[FileLayoutManager]") {
	billiards::unql::UnqlStore<TestClass> manager{"testing.db"};
	std::function<bool(const billiards::unql::Record<TestClass>& r)> printer = [](const billiards::unql::Record<TestClass>& r){
		std::cout << "Entry:\n" << billiards::json::pretty_dump(r) << std::endl;
		return true;
	};
	std::cout << "Initial listing" << std::endl;
	std::cout << "=====================================" << std::endl;
	manager.list(printer);
	std::cout << "=====================================" << std::endl;

	REQUIRE(manager.clear());

	std::cout << "After clearing" << std::endl;
	std::cout << "=====================================" << std::endl;
	manager.list(printer);
	std::cout << "=====================================" << std::endl;

	boost::uuids::uuid uuid{};

	std::string val;
	{
		std::stringstream ss;
		ss << "Here is a number: " << billiards::utils::now();
		val = ss.str();
	}
	{
		billiards::unql::Record<TestClass> record;
		record.obj.str = val;
		manager.create(record);
		uuid = record.info.uuid;
	}

	for (int i = 0; i < 5; i++) {
		billiards::unql::Record<TestClass> record;
		std::stringstream ss;
		ss << "a test thing " << billiards::utils::now() << " " << i << std::endl;
		record.obj.str = ss.str();
		manager.create(record);
	}

	std::cout << "After inserting" << std::endl;
	std::cout << "=====================================" << std::endl;
	manager.list(printer);
	std::cout << "=====================================" << std::endl;

	{
		billiards::unql::Record<TestClass> record;
		bool found = manager.fetch(uuid, record);
		REQUIRE(found);
		REQUIRE(record.obj.str == val);
	}

	manager.remove(uuid);
	std::cout << "After removing " << uuid << std::endl;
	std::cout << "=====================================" << std::endl;
	manager.list(printer);
	std::cout << "=====================================" << std::endl;

	{
		billiards::unql::Record<TestClass> record;
		REQUIRE(!manager.fetch(uuid, record));
	}
}

