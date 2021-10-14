//
// Created by thallock on 10/12/21.
//

#ifndef IDEA_EXTRACTING_FROM_TEST_H
#define IDEA_EXTRACTING_FROM_TEST_H

#include "algebra/parsing.h"
#include "algebra/poly_divide.h"
#include "algebra/VariableNames.h"
#include "algebra/VarietyMatrix.h"
#include "algebra/Variety.h"

namespace billiards::shots::math {

	class LocationsSystem {
	public:
		algebra::poly::VariableNames vars;
		algebra::poly::Variety variety;

		double src_x; double src_y;
		int begin_index;

		inline
		friend std::ostream& operator<<(std::ostream& os, const LocationsSystem& s) {
			std::cout << "Variables (" << s.vars.dim() << "): " << s.vars << std::endl;
			for (int i = 0; i < (int) s.vars.names.size(); i++) {
				if (i != 0) {
					std::cout << ", ";
				}
				std::cout << s.vars.names[i];
			}
			std::cout << std::endl;
			for (int i=0; i < (int) s.variety.polynomials.size(); i++) {
				std::cout << "poly_" << i << " = " << s.vars(s.variety.polynomials[i]) << std::endl;
			}
			return os;
		}

		[[nodiscard]] static inline
		std::string get_var_name(const std::string& type, const std::string& comp, int index) {
			std::stringstream ss;
			ss << type << "_" << index << "_" << comp;
			return ss.str();
		}

		[[nodiscard]] inline
		algebra::poly::Polynomial get_var(const std::string& type, const std::string& comp, int index) const {
			return vars(get_var_name(type, comp, index));
		}
		inline
		void register_var(const std::string& type, const std::string& comp, int index) {
			vars.register_var(get_var_name(type, comp, index));
		}

		[[nodiscard]] inline
		algebra::poly::Polynomial get_exiting_var(const std::string& comp, int step_index) const {
			if (step_index == begin_index - 1) {
				if (comp == "x") {
					return algebra::poly::simplify::constant(vars.dim(), src_x);
				} else if (comp == "y") {
					return algebra::poly::simplify::constant(vars.dim(), src_y);
				} else {
					throw std::runtime_error{"Unknown component"};
				}
			}
			return get_var("s", comp, step_index);
		}

		inline
		void register_exiting_var(const std::string& comp, int step_index) {
			if (step_index == begin_index - 1) {
				return;
			}
			return register_var("s", comp, step_index);
		}

		[[nodiscard]] inline
		algebra::poly::Polynomial get_target_var(const std::string& comp, int step_index) const {
			return get_var("t", comp, step_index);
		}

		inline
		void register_target_var(const std::string& comp, int step_index) {
			return register_var("t", comp, step_index);
		}
	};

	void register_rolling_glance_vars(LocationsSystem& system, int step_index) {
		system.register_exiting_var("x", step_index);
		system.register_exiting_var("y", step_index);
		system.register_target_var("x", step_index);
		system.register_target_var("y", step_index);
		system.register_var("d", "x", step_index);
		system.register_var("d", "y", step_index);
		system.register_var("u", "1", step_index);
		system.register_var("u", "2", step_index);
	}

	void add_glance(
		LocationsSystem& system,
		const double r1, const double r2,
		const double obj_x, const double obj_y,
		int step_index
	) {
		const double alpha = 2.0 / 7;

		auto target_x = system.get_target_var("x", step_index);
		auto target_y = system.get_target_var("y", step_index);
		auto x = system.get_exiting_var("x", step_index);
		auto y = system.get_exiting_var("y", step_index);
		auto ax = system.get_exiting_var("x", step_index - 1);
		auto ay = system.get_exiting_var("y", step_index - 1);
		auto dest_x = system.get_target_var("x", step_index + 1);
		auto dest_y = system.get_target_var("y", step_index + 1);
		auto dx = system.get_var("d", "x", step_index);
		auto dy = system.get_var("d", "y", step_index);
		auto s1 = system.get_var("u", "1", step_index);
		auto s2 = system.get_var("u", "2", step_index);

		// TODO: other assignment as well...
		auto tan_dir_x = (y - obj_y);
		auto tan_dir_y = -(x - obj_x);

		auto tx = x + tan_dir_x;
		auto ty = y + tan_dir_y;

		auto aim_dir_x = x - ax;
		auto aim_dir_y = y - ay;

		auto aim_x = x + s1 * aim_dir_x;
		auto aim_y = y + s1 * aim_dir_y;

		auto destination_1 = x + s2 * (dx - x) - dest_x;
		auto destination_2 = y + s2 * (dy - y) - dest_y;
		auto src_is_exit_1 = target_x - x;
		auto src_is_exit_2 = target_y - y;

		auto glance_1 = aim_x * alpha + tx * (1 - alpha) - dx;
		auto glance_2 = aim_y * alpha + ty * (1 - alpha) - dy;
		auto orth_req_1 = (aim_x - x) * (aim_x - tx) + (aim_y - y) * (aim_y - ty);
		auto radius_1 = (x - obj_x).pow(2) + (y - obj_y).pow(2) - std::pow(r1 + r2, 2);

		system.variety.polynomials.push_back(src_is_exit_1.canon());
		system.variety.polynomials.push_back(src_is_exit_2.canon());
		system.variety.polynomials.push_back(destination_1.canon());
		system.variety.polynomials.push_back(destination_2.canon());
		system.variety.polynomials.push_back(glance_1.canon());
		system.variety.polynomials.push_back(glance_2.canon());
		system.variety.polynomials.push_back(orth_req_1.canon());
		system.variety.polynomials.push_back(radius_1.canon());
	}

	void add_pocket_vars(LocationsSystem& system, int step_index) {
		system.register_target_var("x", step_index);
		system.register_target_var("y", step_index);
	}

	void add_pocket(
		LocationsSystem& system,
		const double r,
		const double px, const double py,
		int step_index
	) {
		auto src_x = system.get_exiting_var("x", step_index - 1);
		auto src_y = system.get_exiting_var("y", step_index - 1);

		auto dest_x = system.get_target_var("x", step_index);
		auto dest_y = system.get_target_var("y", step_index);

		auto orth_req = (dest_x - src_x) * (dest_x - px) + (dest_y - src_y) * (dest_y - py);
		auto radius_req = (dest_x - px).pow(2) + (dest_y - py).pow(2) - std::pow(r, 2);

		system.variety.polynomials.push_back(orth_req.canon());
		system.variety.polynomials.push_back(radius_req.canon());
	}
}

#endif //IDEA_EXTRACTING_FROM_TEST_H
