#include <stdio.h>

#include <functional>
#include <iostream>
#include <cmath>
#include <list>


#define TOL 1e-15
#define LARGER_TOL 1e-8



/////////////////////////////////////////////////////////////////////////////
/// TODO: Solving cubics the hard way...
/////////////////////////////////////////////////////////////////////////////
[[nodiscard]] inline
double eval(const double c0, const double c1, const double c2, const double x) {
	return c0 + c1 * x + c2 * x * x + x * x * x;
}

[[nodiscard]] inline
int get_sign(double value) {
	if (value < -TOL) {
		return -1;
	} else if (value > TOL) {
		return 1;
	} else {
		return 0;
	}
}

[[nodiscard]] inline
bool opposite_signums(const int sign, const int o_sign) {
	return o_sign != 0 && o_sign != sign;
}

inline
void get_sol_bounds(const double c0, const double c1, const double c2, double *lb, double *ub) {
//	int sign = get_sign(eval(c0, c1, c2, 0));
	int sign = get_sign(c0);
	if (sign == 0) {
		throw std::runtime_error{"Simpler case"};
	}
	double t = 1;
	while (t < 1e300) {
		if (opposite_signums(sign, get_sign(eval(c0, c1, c2, t)))) {
			(*lb) = 0;
			(*ub) = t;
			return;
		}
		if (opposite_signums(sign, get_sign(eval(c0, c1, c2, -t)))) {
			(*lb) = 0;
			(*ub) = -t;
			return;
		}
		t *= 2;
	}
	throw std::runtime_error{"This should not happen"};
}

inline
void cubic_binary_search(
	const double c0, const double c1, const double c2,
	const std::function<void(double)>& receiver
) {
	if (std::abs(c0) < TOL) {
		throw std::runtime_error{"constant is zero"};
	} else if (std::abs(c2) < TOL) {
		throw std::runtime_error{"not a cubic"};
	}
	double lb, ub;
	get_sol_bounds(c0, c1, c2, &lb, &ub);
	const int lb_signum = get_sign(eval(c0, c1, c2, lb));
	const int ub_signum = get_sign(eval(c0, c1, c2, ub));
	while (true) {
		const double mid = (lb + ub) / 2.0;
		int signum = get_sign(eval(c0, c1, c2, mid));
		if (signum == 0 || ub - lb < TOL) {
			receiver(mid);
			return;
		}
		if (signum == lb_signum) {
			lb = mid;
			continue;
		}
		if (signum == ub_signum) {
			ub = mid;
			continue;
		}
		throw std::runtime_error{"This should not happen"};
	}
}

[[nodiscard]] inline
double polynomial_divide_13(
	const double c0, const double c1, const double c2, const double c3,
	const double sol,
	double *o0, double* o1, double* o2
) {
	(*o2) = c3;
	(*o1) = c2 + c3 * sol;
	(*o0) = c1 + (c2 + c3 * sol) * sol;
#if 0
	std::cout << "======================================================" << std::endl;
	std::cout << *o0 << " + " << *o1 << " * x + " << *o2 << " * x^2" << std::endl;
	std::cout << "*" << std::endl;
	std::cout << -sol << " + x" << std::endl;
	std::cout << "=" << std::endl;
	std::cout << c0 << " + " << c1 << " * x + " << c2 << " * x^2 + " << c3 << " * x^3" << std::endl;
	std::cout << "======================================================" << std::endl;
#endif
	return c0 + (c1 + (c2 + c3 * sol) * sol) * sol;
	// Solution:
	// 		a3 * x^2 + (a2 + a3 * sol) * x + (a1 + (a2 + a3 * sol) * sol)
	// Multiplication:
	// 		a3 * x^3 + (a2 + a3 * sol) * x^2 + (a1 + (a2 + a3 * sol) * sol) * x
	// 	-	                 a3 * sol  * x^2 +       (a2 + a3 * sol) * sol  * x + (a1 + (a2 + a3 * sol) * sol) * sol
	// --------------------------------------------------------------------------------------------------------------
	// 		a3 * x^3				a2 * x^2       			 		  + a1  * x
}




///////////////////////////////////////////////////////////////////////////////////////////
// TODO:
// FROM: https://www.quora.com/How-can-I-solve-a-cubic-equation-using-C
void solve_cubic_from_website(double a, double b, double c, double d, const std::function<void(double)>& receiver) {
	const double e = 2.7182818284590;
	const double f = ((3 * c / a) - (b * b / (a * a))) / 3; // ** bracketed (a*a)!
	const double g = ((2 * b * b * b / (a * a * a)) - (9 * b * c / (a * a)) + (27 * d / a)) / 27; // ** brackets!
	const double h = (g * g / 4) + (f * f * f / 27);
	const double i = std::sqrt(((g * g / 4) - h));
	const double j = std::exp(std::log10(i) / std::log10(e) / 3);
	const double k = std::acos((-1) * (g / (2 * i)));
	const double l = j * (-1);
	const double m = std::cos(k / 3);
	const double n = std::sqrt(3) * std::sin(k / 3);
	const double p = (b / 3 * a) * (-1);
	const double r = (-1) * (g / 2) + std::sqrt(h);
	const double s = std::exp(log10(r) / std::log10(e) / 3);
	const double t = (-1) * (g / 2) - std::sqrt(h);
	const double u = std::exp(log10(t) / std::log10(e) / 3);

	bool w_assigned = false;
	int w;
	if (h > 0) {
		w = 1;
		w_assigned = true;
	}
	if (h <= 0) {
		w = 3;
		w_assigned = true;
	}
	if (std::abs(f) < TOL && std::abs(g) < TOL && std::abs(h) < TOL) {
		w = 2;
		w_assigned = true;
	}
	if (!w_assigned) {
		// TODO: Not sure if this can happen...
		return;
	}
	switch (w) {
		case 1: {
			const double x1 = (s + u) - (b / 3 * a);
			const double x2 = (-1) * (s + u) / 2 - (b / 3 * a);
			const double x3 = (s - u) * sqrt(3) / 2;
			receiver(x1);
			// sol 2: x2 + i * x3
			// sol 3: x2 - i * x3
			break;
		}
		case 2: {
			const double x = std::exp(std::log10(d / a) / std::log10(e) / 3) * (-1);
//			printf("\n There is a line:\n%lf", x1);
			break;
		}
		case 3: {
			receiver(2 * j * cos(k / 3) - (b / 3 * a));
			receiver(l * (m + n) + p);
			receiver(l * (m - n) + p);
			break;
		}
	}
}



/////////////////////////////////////////////////////////////////////////////
/// TODO: Polynomials...
/////////////////////////////////////////////////////////////////////////////

// Solve c == 0
inline
void solve_0(
	const double c,
	const double default_value,
	const std::function<void(const double)>& receiver
) {
	if (std::abs(c) < TOL) {
		receiver(default_value);
	}
}

// Solve c + x == 0
inline
void solve_1n(
	const double c,
	const std::function<void(const double)>& receiver
) {
	receiver(-c);
}

// Solve c0 + c1 * x == 0
inline
void solve_1(
	const double c0, const double c1,
	const double default_value,
	const std::function<void(const double)>& receiver
) {
	if (std::abs(c1) < TOL) {
		solve_0(c0, default_value, receiver);
		return;
	}
	solve_1n(c0 / c1, receiver);
}

// Solve c0 + c1 * x + x^2 == 0
inline
void solve_2n(
	const double c0, const double c1,
	const std::function<void(const double)>& receiver
) {
	const double discriminant = c1 * c1 - 4 * c0;
	if (discriminant < 0) {
		return;
	}
	if (discriminant < TOL) {
		receiver(-c1 / 2);
		return;
	}
	const double d = std::sqrt(discriminant);
	receiver((-c1 + d) / 2);
	receiver((-c1 - d) / 2);
}

// Solve c0 + c1 * x + c2 * x^2 == 0
inline
void solve_2(
	const double c0, const double c1, const double c2,
	const double default_value,
	const std::function<void(const double)>& receiver
) {
	if (std::abs(c2) < TOL) {
		solve_1(c0, c1, default_value, receiver);
		return;
	}
	solve_2n(c0 / c2, c1 / c2, receiver);
}


// Solve c0 + c1 * x + c2 * x^2 + x^3 == 0, the hard way
inline
void solve_3n_binary_search(
	const double c0, const double c1, const double c2,
	const double default_value,
	const std::function<void(const double)>& receiver
) {
	if (std::abs(c0) < TOL) {
		receiver(0);
		solve_2n(c1, c2, receiver);
		return;
	}
	cubic_binary_search(
		c0, c1, c2,
		[c0, c1, c2, default_value, &receiver] (double sol) {
			receiver(sol);

			double a0, a1, a2;
			const double r = polynomial_divide_13(
				c0, c1, c2, 1.0,
				sol,
				&a0, &a1, &a2
			);
			std::cout << "Remainder polynomial: " << a0 << ", " << a1 << ", " << a2 << std::endl;
			std::cout << "Remainder: " << r << std::endl;
			solve_2(a0, a1, a2, default_value, receiver);
		}
	);
}




#if 0
inline
void solve_3n(
	const double c, const double b, const double a,
	const double default_value,
	const std::function<void(const double)>& receiver
) {
	if (std::abs(c) < TOL) {
		receiver(0);
		solve_2n(b, a, receiver);
		return;
	}

	const double p = (3 * b - a * a) / 3;
	const double q = (2 * a * a * a - 9 * a * b + 27 * c) / 27;
	const double r1 = q * q + p * p * p / 27;
	const double r2 = q * q / 4 + p * p * p / 27;
	if (r1 < 0 || r2 < 0) {
//		solve_3n_binary_search(c, b, a, default_value, receiver);
		return;
	}
	const double sol = std::pow(-q / 2 + std::sqrt(r1), 1/3.0) + std::pow(-q / 2 - std::sqrt(r2), 1/3.0) - a / 3;
	receiver(sol);

	double a0, a1, a2;
	const double r = polynomial_divide_13(
		c, b, a, 1.0,
		sol,
		&a0, &a1, &a2
	);
	solve_2(a0, a1, a2, default_value, receiver);
}
#endif


class Complex {
public:
	double a;
	double b;

	Complex(double a, double b) : a{a}, b{b} {}
	Complex(double a) : a{a}, b{0} {}

#define FROM_EXP(r, theta) Complex{r * std::cos(theta), r * std::sin(theta)}

	[[nodiscard]] inline
	double r() const {
		return std::sqrt(a * a + b * b);
	}
	[[nodiscard]] inline
	double theta() const {
		return std::atan2(b, a);
	}
	[[nodiscard]] inline
	bool is_real() const {
		return std::abs(b) < LARGER_TOL;
	}
	[[nodiscard]] inline
	double real() const {
		return a;
	}
	[[nodiscard]] inline
	Complex operator*(const Complex& o) const {
		return Complex{a * o.a - b * o.b, a * o.b + b * o.a};
	}
	[[nodiscard]] inline
	Complex operator+(const Complex& o) const {
		return Complex{a + o.a, b + o.b};
	}
	[[nodiscard]] inline
	Complex operator-(const Complex& o) const {
		return Complex{a - o.a, b - o.b};
	}
	[[nodiscard]] inline
	Complex operator/(const Complex& o) const {
		return FROM_EXP(r() / o.r(), theta() - o.theta());
	}
	[[nodiscard]] inline
	Complex operator*(const double d) const {
		return Complex{a * d, b * d};
	}
	[[nodiscard]] inline
	Complex operator+(const double d) const {
		return Complex{a + d, b};
	}
	[[nodiscard]] inline
	Complex pow(const double exp) const {
		return FROM_EXP(std::pow(r(), exp), theta() * exp);
	}
	friend std::ostream& operator<<(std::ostream& os, const Complex& c) {
		return os << c.a << " + " << c.b << "i";
	}
};

// Solve -Q + P * x + x^3 = 0
inline
void solve_3_simp(
	const double P, const double Q,
	const std::function<void(const double)>& receiver
) {
	// https://www.shsu.edu/kws006/professional/Concepts_files/SolvingCubics.pdf
	const double delta = std::pow(P / 3, 3) + std::pow(Q / 2, 2);
	const Complex b = Complex{delta}.pow(0.5) + (-Q / 2);
	const Complex beta = b.pow(1 / 3.0);
	const Complex alpha = Complex{P} / (beta * 3);
	const Complex omega = Complex{-0.5, std::sqrt(3) / 2};
	for (const Complex& sol : std::array<Complex, 3>{
		alpha - beta,
		alpha * omega * omega - beta * omega,
		alpha * omega - beta * omega * omega
	}) {
		if (sol.is_real()) {
			receiver(sol.real());
		}
	}
}

/*
c0, c1, c2, x = var('c0 c1 c2 x')
y = x - c2 / 3
expand(c0 + c1 * y + c2 * y^2 + y^3).collect(x)
2/27*c2^3 + x^3 - 1/3*c1*c2 - 1/3*(c2^2 - 3*c1)*x + c0
*/
inline
void solve_3n(
	const double c0, const double c1, const double c2,
	const double default_value,
	const std::function<void(const double)>& receiver
) {
	const double P = -1/3.0*(c2*c2 - 3*c1);
	const double Q = -(2/27.0*c2*c2*c2 - 1/3.0*c1*c2 + c0);

	std::cout << "===========" << std::endl;
	std::cout << "ayo" << std::endl;
	const double vx = 24;
	const double vy = vx - c2 / 3;
	std::cout << c0 + c1 * vy + c2 * vy * vy + vy * vy * vy << std::endl;
	std::cout << (2/27.0*c2*c2*c2 + vx*vx*vx - 1/3.0*c1*c2 - 1/3.0*(c2*c2 - 3*c1)*vx + c0) << std::endl;
	std::cout << (2/27.0*c2*c2*c2 + -1/3.0*c1*c2 + c0) << ", " << -Q << std::endl;
	std::cout <<  - 1/3.0*(c2*c2 - 3*c1) << ", " << P << std::endl;
	std::cout << vx * vx * vx + P * vx - Q << std::endl;
	std::cout << "===========" << std::endl;

//	const double P = -c2 * c2 / 3;
//	const double Q = -(c0 + 2 * c2 * c2 * c2 / 27 - c1 * c2 / 3);

	solve_3_simp(
		-1/3.0*(c2*c2 - 3*c1), -(2/27.0*c2*c2*c2 - 1/3.0*c1*c2 + c0),
		[P, Q, c0, c1, c2, &receiver](double x) {
			std::cout << "Solves it:" << std::endl;
			std::cout << (x * x * x + P * x - Q) << std::endl;
			const double y = x - c2 / 3;
			std::cout << c0 + c1 * y + c2 * y * y + y * y * y << std::endl;
			std::cout << "---" << std::endl;
			receiver(x - c2 / 3);
		});
}





// Solve c0 + c1 * x + c2 * x^2 + c3 * x^3 == 0
inline
void solve_3(
	const double c0, const double c1, const double c2, const double c3,
	const double default_value,
	const std::function<void(const double)>& receiver
) {
	if (std::abs(c3) < TOL) {
		solve_2(c0, c1, c2, default_value, receiver);
		return;
	}
	solve_3n(c0 / c3, c1 / c3, c2 / c3, default_value, receiver);
}

// Solve:
// a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + x ** 4 == 0
void solve_4n(
	const double a0, const double a1, const double a2, const double a3,
	const double default_val,
	const std::function<void(const double)>& receiver
) {
	std::cout << "Calling it" << std::endl;
	// https://mathworld.wolfram.com/QuarticEquation.html

	std::list<double> found_solutions;
	const std::function<void(const double)> check_receiver = [a0, a1, a2, a3, &receiver, &found_solutions](const double s) {
		if (std::abs(a0 + a1 * s + a2 * s * s + a3 * s * s * s + s * s * s * s) > LARGER_TOL) {
			return;
		}
		for (const double f : found_solutions) {
			if (std::abs(s - f) < LARGER_TOL) {
				return;
			}
		}
		found_solutions.push_back(s);
		receiver(s);
	};
	solve_3n(
		4 * a2 * a0 - a1 * a1 - a3 * a3 * a0, a1 * a3 - 4 * a0, -a2,
		default_val,
		[a0, a1, a2, a3, &check_receiver](const double y1) {
			const double r1 = a3 * a3 - 4 * a2 + 4 * y1;
			const double r2 = y1 * y1 - 4 * a0;
			if (r1 < 0 || r2 < 0) {
				return;
			}

			std::cout << "Solves the cubic: " << std::endl;
			std::cout << (4 * a2 * a0 - a1 * a1 - a3 * a3 * a0 + (a1 * a3 - 4 * a0) * y1 - a2 * y1 * y1 + y1 * y1 * y1) << std::endl;


			solve_2n(0.5 * (y1 + std::sqrt(r2)), 0.5 * (a3 + std::sqrt(r1)), [y1, a0, a1, a2, a3, a3, &check_receiver, r1, r2](double y2) {
				std::cout << "quadratic solve: " << (0.5 * (y1 + std::sqrt(r2)) + 0.5 * (a3 + std::sqrt(r1)) * y2 + y2 * y2) << std::endl;
				std::cout << "quartic sol: " << (
					a0 + a1 * y2 + a2 * y2 * y2 + a3 * y2 * y2 * y2 + y2 * y2 * y2 * y2
					) << std::endl;
				check_receiver(y2);
			});
			solve_2n(0.5 * (y1 - std::sqrt(r2)), 0.5 * (a3 - std::sqrt(r1)), [y1, a0, a1, a2, a3, &check_receiver, r1, r2](double y2) {
				std::cout << "quadratic solve: " << (0.5 * (y1 - std::sqrt(r2)) + 0.5 * (a3 - std::sqrt(r1)) * y2 + y2 * y2) << std::endl;
				std::cout << "quartic sol: " << (
					a0 + a1 * y2 + a2 * y2 * y2 + a3 * y2 * y2 * y2 + y2 * y2 * y2 * y2
				) << std::endl;
				check_receiver(y2);
			});
			solve_2n(0.5 * (y1 - std::sqrt(r2)), 0.5 * (a3 + std::sqrt(r1)), [y1, a0, a1, a2, a3, &check_receiver, r1, r2](double y2) {
				std::cout << "quadratic solve: " << (0.5 * (y1 - std::sqrt(r2)) + 0.5 * (a3 + std::sqrt(r1)) * y2 + y2 * y2) << std::endl;
				std::cout << "quartic sol: " << (
					a0 + a1 * y2 + a2 * y2 * y2 + a3 * y2 * y2 * y2 + y2 * y2 * y2 * y2
				) << std::endl;
				check_receiver(y2);
			});
			solve_2n(0.5 * (y1 + std::sqrt(r2)), 0.5 * (a3 - std::sqrt(r1)), [y1, a0, a1, a2, a3, &check_receiver, r1, r2](double y2) {
				std::cout << "quadratic solve: " << (0.5 * (y1 + std::sqrt(r2)) + 0.5 * (a3 - std::sqrt(r1)) * y2 + y2 * y2) << std::endl;
				std::cout << "quartic sol: " << (
					a0 + a1 * y2 + a2 * y2 * y2 + a3 * y2 * y2 * y2 + y2 * y2 * y2 * y2
				) << std::endl;
				check_receiver(y2);
			});
		}
	);
}

// Solve:
// a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + a4 * x ** 4 == 0
void solve_4(
	const double a0, const double a1, const double a2, const double a3, const double a4,
	const double default_val,
	const std::function<void(const double)>& sol
) {
	if (std::abs(a4) < TOL) {
		solve_3(a0, a1, a2, a3, default_val, sol);
		return;
	}
	solve_4n(a0 / a4, a1 / a4, a2 / a4, a3 / a4, default_val, sol);
}


int main() {
	if (0) {
		Complex c{3.4, -2.9};
		Complex r = c.pow(1 / 3.0);
		std::cout << c << std::endl;
		std::cout << r << std::endl;
		std::cout << (r * r * r) << std::endl;
/*
x = var('x')
expand((x - 3) * (x + 2) * (x-7) * (x + 4))
expand((x - 3) * (x + 2) * (x-7))

 */
		const double c4 = 1;
		const double c3 = -4;
		const double c2 = -31;
		const double c1 = 46;
		const double c0 = 168;

		{
			for (const double x: std::array<double, 4>{3, -2, 7, -4}) {
				std::cout << (c0 + c1 * x + c2 * x * x + c3 * x * x * x + c4 * x * x * x * x) << std::endl;
			}
			std::cout << "-" << std::endl;
			for (const double x: std::array<double, 3>{3, -2, 7}) {
				std::cout << (42 + x - 8 * x * x + x * x * x) << std::endl;
			}
			std::cout << "-" << std::endl;
		}

		const double P = -6;
		const double Q = 4;
		solve_3_simp(P, Q, [P, Q](double x) {
			std::cout << "Simplified cubic solution: " << x << std::endl;
			std::cout << "value: " << (x * x * x + P * x - Q) << std::endl;
		});

		solve_2(
			12, -7, 1.0,
			0,
			[](double x) {
				std::cout << "Quadratic solution: " << x << std::endl;
				std::cout << "p(x) = " << (12 - 7 * x + x * x) << std::endl;
			});

		solve_3(
			42, 1.0, -8, 1.0,
			0,
			[](double x) {
				std::cout << "Cubic solution: " << x << std::endl;
				std::cout << "p(x) = " << (42 + x - 8 * x * x + x * x * x) << std::endl;
			});

#if 0
		std::cout << "Sovling:" << std::endl;
		std::cout << "p(x) = " << c3 << "*x^3 + " << c2 << "*x^2 + " << c1 << "*x + " << c0 << " == 0" << std::endl;
		solve_cubic_from_website(
			c3, c2, c1, c0,
			[c0, c1, c2, c3](double x) {
				std::cout << "Website solution x = " << x << std::endl;
				std::cout << "p(x) = " << (c0 + c1 * x + c2 * x * x + c3 * x * x * x) << std::endl;
			}
		);
#endif
		solve_4(
			c0, c1, c2, c3, c4,
			0,
			[c0, c1, c2, c3, c4](double x) {
				std::cout << "Quartic solution: " << x << std::endl;
				std::cout << "p(x) = " << (c0 + c1 * x + c2 * x * x + c3 * x * x * x + c4 * x * x * x * x) << std::endl;
			});
	}

	{
		const double c0 = -0.034133333333333335;
		const double c1 = 0.17066666666666666;
		const double c2 = 0;
		const double c3 = -0.53333333333333333;
		// -0.034133333333333335 + 0.17066666666666666 * x + -0.53333333333333333 * x^2 + x^4

		solve_4n(
			c0, c1, c2, c3,
			0,
			[c0, c1, c2, c3](double x) {
				std::cout << "Quartic solution: " << x << std::endl;
				std::cout << "p(x) = " << (c0 + c1 * x + c2 * x * x + c3 * x * x * x + x * x * x * x) << std::endl;
			});
	}

	return 0;
}











#if 0


// q^2 = 4 * (u - p) * (u^2 / 4 - r)
// q^2 = 4 * (u - p) * u^2 / 4 - 4 * (u - p) * r
// q^2 = 4 * u * u^2 / 4 - 4 * p * u^2 / 4 - 4 * r * u + 4 * r * p
// q^2 = u^3 - 4 * p * u^2 - 4 * r * u + 4 * r * p
// 0 = u^3 - (q^2 + 4 * p) * u^2 - 4 * r * u + 4 * r * p


// Solve:
// a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + x ** 4 == 0
void solve_4n(
	const double c0, const double c1, const double c2, const double c3,
	const double default_val,
	const std::function<void(const double)>& receiver
) {
	// https://www.maa.org/sites/default/files/pdf/upload_library/22/Ford/auckly29.pdf
	// d == c0
	// c == c1
	// b == c2
	// a == c3

	const double p = c2 - 3 * c3 * c3 / 8;
	const double q = c1 - c3 * c2 / 2 + c3 * c3 * c3 / 8;
	const double r = c0 - c3 * c1 / 4 + c3 * c3 * c2 / 16 - 3 * c3 * c3 * c3 * c3 / 256;
	solve_3n(
//		q * q, p * p - 4 * r, -2 * p,
		-q * q, p * p - 4 * r, 2 * p,
//		-q * q, p * p - 4 * r, -2 * p,
//		q * q, p * p - 4 * r, -2 * p,
		default_val,
		[p, q, c3, receiver](const double lam){
			if (lam < -TOL) {
				// TODO
				return;
			}
			if (lam < TOL) {
				// TODO
				return;
			}
			const double slam = std::sqrt(lam);

			const double prad = lam - 2 * (p + lam + q / slam);
			if (prad > 0) {
				receiver((-slam + std::sqrt(prad)) / 2 - c3 / 4);
				receiver((-slam - std::sqrt(prad)) / 2 - c3 / 4);
			}
			const double nrad = lam - 2 * (p + lam - q / slam);
			if (nrad > 0) {
				receiver((+slam + std::sqrt(nrad)) / 2 - c3 / 4);
				receiver((+slam - std::sqrt(nrad)) / 2 - c3 / 4);
			}
		}
	);
}

#endif