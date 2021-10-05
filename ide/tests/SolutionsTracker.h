//
// Created by thallock on 10/4/21.
//

#ifndef IDEA_SOLUTIONSTRACKER_H
#define IDEA_SOLUTIONSTRACKER_H

namespace billiards::test {

	class SolutionsTracker {
	public:
		std::shared_ptr<std::vector<double>> sols;

		SolutionsTracker(const SolutionsTracker& tracker) = default;
		SolutionsTracker() : sols{std::make_shared<std::vector<double>>()} {}
		~SolutionsTracker() = default;

		inline
		void operator()(const double value) const {
			sols->push_back(value);
		}

		[[nodiscard]] inline
		int count() const {
			return (int) sols->size();
		}

		[[nodiscard]] inline
		double operator[](const int idx) const {
			return sols->at(idx);
		}

		[[nodiscard]] inline
		bool contains(const double value) const {
			for (const double sol : *sols) {
				if (std::abs(sol - value) < 1e-6) {
					return true;
				}
			}
			return false;
		}
	};


	class SolutionsTracker2d {
	public:
		std::shared_ptr<std::vector<std::pair<double, double>>> sols;

		SolutionsTracker2d(const SolutionsTracker2d& tracker) : sols{tracker.sols} {}
		SolutionsTracker2d() : sols{std::make_shared<std::vector<std::pair<double, double>>>()} {}
		~SolutionsTracker2d() = default;

		inline
		void operator()(const double x, const double y) const {
			sols->emplace_back(x, y);
		}

		[[nodiscard]] inline
		int count() const {
			return (int) sols->size();
		}

		[[nodiscard]] inline
		std::pair<double, double> operator[](const int idx) const {
			return sols->at(idx);
		}

		[[nodiscard]] inline
		bool contains(const std::pair<double, double>& value) const {
			for (const std::pair<double, double>& sol : *sols) {
				if (std::abs(sol.first - value.first) < LARGER_TOL
					&& std::abs(sol.second - value.second) < LARGER_TOL
				) {
					return true;
				}
			}
			return false;
		}
	};
}


#endif //IDEA_SOLUTIONSTRACKER_H
