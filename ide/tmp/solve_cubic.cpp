#include <stdio.h>

#include <functional>
#include <iostream>
#include <cmath>
#include <list>


#define TOL 1e-15
#define LARGER_TOL 1e-7



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
	std::cout << "P: " << P << std::endl;
	std::cout << "Q: " << Q << std::endl;

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

// Solves:
// 		c0 + c1 * x + c2 * x * x + x * x * x
inline
void solve_3n(
	const double c0, const double c1, const double c2,
	const std::function<void(const double)>& receiver
) {
	if (std::abs(c0) < TOL) {
		receiver(0);
		solve_2n(c1, c2, receiver);
	}

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
	std::cout << "c0: " << c0 << std::endl;
	std::cout << "c1: " << c1 << std::endl;
	std::cout << "c2: " << c2 << std::endl;
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
	solve_3n(c0 / c3, c1 / c3, c2 / c3, receiver);
}

// Solve:
// a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3 + x ** 4 == 0
void solve_4n(
	const double a0, const double a1, const double a2, const double a3,
	const std::function<void(const double)>& receiver
) {
	if (std::abs(a0) < TOL) {
		receiver(0);
		solve_3n(a1, a2, a3, receiver);
	}
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
		[a0, a1, a2, a3, &check_receiver](const double y1) {
			const double r1a = a3 * a3 - 4 * a2 + 4 * y1;
			const double r2a = y1 * y1 - 4 * a0;
			const double r1 = std::abs(r1a) < LARGER_TOL ? 0 : r1a;
			const double r2 = std::abs(r2a) < LARGER_TOL ? 0 : r2a;
			std::cout << "r1=" << r1 << ", r2=" << r2 << std::endl;
			if (r1 < 0 || r2 < 0) {
				return;
			}

			std::cout << "Solves the cubic: " << std::endl;
			std::cout << (4 * a2 * a0 - a1 * a1 - a3 * a3 * a0 + (a1 * a3 - 4 * a0) * y1 - a2 * y1 * y1 + y1 * y1 * y1) << std::endl;


			solve_2n(0.5 * (y1 + std::sqrt(r2)), 0.5 * (a3 + std::sqrt(r1)), [y1, a0, a1, a2, a3, &check_receiver, r1, r2](double y2) {
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
	solve_4n(a0 / a4, a1 / a4, a2 / a4, a3 / a4, sol);
}


// Solve:
// a00 + a10 * x + a20 * x ** 2 == 0
void solve_20(
	const double a00, const double a10, const double a20,
	const double default_x, const double default_y,
	const std::function<void(const double x, const double y)>& receiver
) {
	// Could use default_y
	solve_2(a00, a10, a20, default_x, [default_y, &receiver](const double x) { receiver(x, default_y); });
}

// Solve
// a00 + a10 * x + a01 * y + a20 * x ** 2 == 0
void solve_21(
	const double a00, const double a10, const double a01, const double a20,
	const double default_x, const double default_y,
	const std::function<void(const double x, const double y)>& receiver
) {
	if (std::abs(a01) < TOL) {
		solve_20(a00, a10, a20, default_x, default_y, receiver);
		return;
	}
	solve_1n(
		(a00 + a10 * default_x + a20 * default_x * default_x) / a01,
		[default_x, receiver](const double y) {
			receiver(default_x, y);
		});
}

// Solve
// a00 + a10 * x + a01 * y + a20 * x ** 2 + a02 * y ** 2 == 0
void solve_22_simple(
	const double a00, const double a10, const double a01, const double a20, const double a02,
	const double default_x, const double default_y,
	const std::function<void(const double x, const double y)>& receiver
) {
	if (std::abs(a02) < TOL) {
		solve_21(a00, a10, a01, a20, default_x, default_y, receiver);
		return;
	}
	solve_2n(
		(a00 + a10 * default_x + a20 * default_x * default_x) / a02,
		a01 / a02,
		[default_x, &receiver](const double y) {
			receiver(default_x, y);
		}
	);
}

// Solve:
// a00 + a10 * x + a01 * y + a20 * x ** 2 + a11 * x * y + a02 * y ** 2 == 0
void solve_22(
	const double a00, const double a10, const double a01, const double a20, const double a11, const double a02,
	const double default_x, const double default_y,
	const std::function<void(const double x, const double y)>& receiver
) {
	if (std::abs(a11) < TOL) {
		solve_22_simple(a00, a10, a01, a20, a02, default_x, default_y, receiver);
		return;
	}

	// try default_x...
	{
		const double b0 = a00 + a10 * default_x + a20 * default_x * default_x;
		const double b1 = a01 + a11 * default_x;
		const double b2 = a02;

		if (std::abs(b2) > TOL) {
			solve_2n(
				b0 / b2, b1 / b2,
				[default_x, &receiver](const double y) {
					receiver(default_x, y);
				}
			);
			return;
		}
		if (std::abs(b1) > TOL) {
			solve_1n(
				b0 / b1,
				[default_x, &receiver](const double y) {
					receiver(default_x, y);
				}
			);
		}
	}

	// try default_y...
	// a00 + a10 * x + a01 * y + a20 * x ** 2 + a11 * x * y + a02 * y ** 2 == 0
	{
		const double b0 = a00 + a01 * default_y + a02 * default_y * default_y;
		const double b1 = a10 + a11 * default_y;
		const double b2 = a20;

		if (std::abs(b2) > TOL) {
			solve_2n(
				b0 / b2, b1 / b2,
				[default_y, &receiver](const double x) {
					receiver(x, default_y);
				}
			);
			return;
		}
		if (std::abs(b1) > TOL) {
			solve_1n(
				b0 / b1,
				[default_y, &receiver](const double x) {
					receiver(x, default_y);
				}
			);
		}
	}
	// a02 == 0
	// a20 == 0
	// a11 =/= 0
	// a10 + a11 * default_y == 0
	// a01 + a11 * default_x == 0

	// a00 - a11 * default_y * x - a11 * default_x * y + a11 * x * y == 0
	// a00 + a11 * (-default_y * x - default_x * y + x * y) == 0
	throw std::runtime_error{"Implement this case: Infinitely many solutions exist, but default_x, default_y does not work"};
}

// Solve:
// a00 == 0
// b00 + b10 * x + b01 * y + b20 * x ** 2 + b11 * x * y + b02 * y ** 2 == 0
void solve_00_22(
	const double a00,
	const double b00, const double b10, const double b01, const double b20, const double b11, const double b02,
	const double default_x, const double default_y,
	const std::function<void(const double x, const double y)>& receiver
) {
	if (std::abs(a00) < TOL) {
		solve_22(b00, b10, b01, b20, b11, b02, default_x, default_y, receiver);
		return;
	}
}

// Substitute x = x_val into
// a00 + a10 * x + a01 * y + a20 * x ** 2 + a11 * x * y + a02 * y ** 2 == 0
void subs_22_x(
	const double x_val,
	const double a00, const double a10, const double a01, const double a20, const double a11, const double a02,
	const double default_y,
	const std::function<void(const double x, const double y)>& receiver
) {
	solve_2(
		a00 + a10 * x_val + a20 * x_val * x_val, a01 + a11 * x_val, a02,
		default_y,
		[x_val, &receiver](const double y) { receiver(x_val, y); }
	);
}
// Substitute y = y_val into
// a00 + a10 * x + a01 * y + a20 * x ** 2 + a11 * x * y + a02 * y ** 2 == 0
void subs_22_y(
	const double y_val,
	const double a00, const double a10, const double a01, const double a20, const double a11, const double a02,
	const double default_x,
	const std::function<void(const double x, const double y)>& receiver
) {
	solve_2(
		a00 + a01 * y_val + a02 * y_val * y_val, a10 + a11 * y_val, a20,
		default_x,
		[y_val, &receiver](const double x) { receiver(x, y_val); }
	);
}

// Solve:
// a00 + a10 * x == 0
// b00 + b10 * x + b01 * y + b20 * x ** 2 + b11 * x * y + b02 * y ** 2 == 0
void solve_01_22(
	const double a00, const double a10,
	const double b00, const double b10, const double b01, const double b20, const double b11, const double b02,
	const double default_x, const double default_y,
	const std::function<void(const double x, const double y)>& receiver
) {
	if (std::abs(a10) < TOL) {
		solve_00_22(a00, b00, b10, b01, b20, b11, b02, default_x, default_y, receiver);
		return;
	}

	solve_1n(
		a00 / a10,
		[b00, b01, b10, b20, b11, b02, default_y, &receiver](const double x) {
			subs_22_x(x, b00, b10, b01, b20, b11, b02, default_y, receiver);
		}
	);
}

// Solve:
// a00 + a10 * x + a20 * x ** 2 == 0
// b00 + b10 * x + b01 * y + b20 * x ** 2 + b11 * x * y + b02 * y ** 2 == 0
void solve_02_22(
	const double a00, const double a10, const double a20,
	const double b00, const double b10, const double b01, const double b20, const double b11, const double b02,
	double default_x, double default_y,
	const std::function<void(const double, const double)>& receiver
) {
	if (std::abs(a20) < TOL) {
		solve_01_22(a00, a10, b00, b10, b01, b20, b11, b02, default_x, default_y, receiver);
		return;
	}

	solve_2n(
		a00 / a20, a10 / a20,
		[b00, b10, b01, b20, b11, b02, default_y, &receiver](const double x) {
			subs_22_x(x, b00, b10, b01, b20, b11, b02, default_y, receiver);
		}
	);
}

// Solve:
// a00 + a10 * x + a01 * y + a20 * x ** 2 == 0
// b00 + b10 * x + b01 * y + b20 * x ** 2 + b11 * x * y + b02 * y ** 2 == 0
/*
 *
x, y = var('x y')
a00, a10, a01, a20 = var('a00 a10 a01 a20')
b00, b10, b01, b20, b11, b02 = var('b00 b10 b01 b20 b11 b02')
eq1 = a00 + a10 * x + a01 * y + a20 * x ** 2
eq2 = b00 + b10 * x + b01 * y + b20 * x ** 2 + b11 * x * y + b02 * y ** 2

y_sol = solve(eq1 == 0, y)[0].right()
poly = eq2.substitute(y=y_sol).expand().collect(x)

print('\t\tconst double y = ' + convert_powers(y_sol) + ';')

for coef, ord in poly.coefficients(x):
	print('\t\tconst double q' + str(ord) + ' = ' + convert_powers(coef) + ';')

 */
void solve_21_22(
	const double a00, const double a10, const double a01, const double a20,
	const double b00, const double b01, const double b10, const double b20, const double b11, const double b02,
	double default_x, double default_y,
	const std::function<void(const double, const double)>& receiver
) {
	if (std::abs(a01) < TOL) {
		solve_02_22(a00, a10, a20, b00, b10, b01, b20, b11, b02, default_x, default_y, receiver);
		return;
	}

	const double q0 = b00 - a00*b01/a01 + std::pow(a00, 2)*b02/std::pow(a01, 2);
	const double q1 = -a10*b01/a01 + 2*a00*a10*b02/std::pow(a01, 2) + b10 - a00*b11/a01;
	const double q2 = -a20*b01/a01 + std::pow(a10, 2)*b02/std::pow(a01, 2) + 2*a00*a20*b02/std::pow(a01, 2) - a10*b11/a01 + b20;
	const double q3 = 2*a10*a20*b02/std::pow(a01, 2) - a20*b11/a01;
	const double q4 = std::pow(a20, 2)*b02/std::pow(a01, 2);
	solve_4(
		q0, q1, q2, q3, q4,
		default_x
		[a20, a10, a00, receiver](const double x) {
			const double y = -(a20*std::pow(x, 2) + a10*x + a00)/a01;
			receiver(x, y);
		}
	);
}


// Solve:
// a0 * x^2 + p0 * x + q0 == 0
// a1 * x^2 + p1 * x + q1 == 0
void solve_20_20(
	const double q0, const double p0, const double a0,
	const double q1, const double p1, const double a1,
	const double default_x,
	const std::function<void(const double)> receiver
) {
	const auto check0_rec = [a0, p0, q0, &receiver](const double x) {
		if (std::abs(a0 * x * x + p0 * x + q0) < LARGER_TOL) {
			receiver(x);
		}
	};
	const auto check1_rec = [a1, p1, q1, &receiver](const double x) {
		if (std::abs(a1 * x * x + p1 * x + q1) < LARGER_TOL) {
			receiver(x);
		}
	};
	if (std::abs(a0) > TOL) {
		solve_2n(q0 / a0, p0 / a0, check1_rec);
		return;
	}
	if (std::abs(a1) > TOL) {
		solve_2n(q1 / a1, p1 / a1, check0_rec);
		return;
	}
	if (std::abs(p0) > TOL) {
		solve_1n(q0 / p0, check1_rec);
	}
	if (std::abs(p1) > TOL) {
		solve_1n(q1 / p1, check0_rec);
	}
	if (std::abs(q0) < TOL && std::abs(q1) < TOL) {
		return;
	}
	receiver(default_x);

// a0 * x^2 + p0 * x + q0 == 0
// a1 * x^2 + p1 * x + q1 == 0
}


/*
 * Sage script:

import re
def convert_powers(expr):
	return re.sub(r'([a-zA-Z0-9]*)\^([0-9]+)', r'std::pow(\g<1>, \g<2>)', str(expr))

# R.<x,y> = QQ[]
x, y = var('x y')
a0, b0, c0, d0, e0, f0 = var('a0 b0 c0 d0 e0 f0')
a1, b1, c1, d1, e1, f1 = var('a1 b1 c1 d1 e1 f1')

eq00 = f0 + d0 * x + e0 * y + a0 * x^2 + c0 * x * y + b0 * y^2
eq01 = f1 + d1 * x + e1 * y + a1 * x^2 + c1 * x * y + b1 * y^2

p0 = c0 * y + d0
q0 = f0 + e0 * y + b0 * y^2
p1 = c1 * y + d1
q1 = f1 + e1 * y + b1 * y^2

eq10 = a0 * x^2 + p0 * x + q0
eq11 = a1 * x^2 + p1 * x + q1

bool(eq00 == eq10) # True
bool(eq01 == eq11) # True

eq2 = (a1 * eq10 - a0 * eq11).expand().collect(x)

for coef, ord in eq2.coefficients(x):
	print('\t\t\tconst double l' + str(ord) + ' = ' + convert_powers(coef) + ';')

# x_sol = eq2.roots(x)[0][0]
x_sol = solve(eq2, x)[0].right()
num = numerator(x_sol)
den = denominator(x_sol)
l0, l1 = [p[0] for p in sorted(eq2.collect(x).coefficients(x), key=lambda p: p[1])]
bool(l0 + l1 * x == eq2)
bool(x_sol == num / den)

# print('\t\t\tconst double num = ' + convert_powers(num) + ';')
# print('\t\t\tconst double den = ' + convert_powers(den) + ';')

eq3 = (den^2 * eq10.substitute(x=x_sol)).simplify_full().expand().collect(y)
# eq3 = (a0 * num^2 + p0 * num * den + q0 * den^2).collect(y)

for coef, ord in eq3.coefficients(y):
	print('\t\tconst double q' + str(ord) + ' = ' + convert_powers(coef) + ';')

 */
//
// Solve:
// 		f0 + d0 * x + e0 * y + a0 * x ** 2 + c0 * x * y + b0 * y ** 2 == 0
// 		f1 + d1 * x + e1 * y + a1 * x ** 2 + c1 * x * y + b1 * y ** 2 == 0
//
void solve_22_22(
	const double f0, const double d0, const double e0, const double a0, const double c0, const double b0,
	const double f1, const double d1, const double e1, const double a1, const double c1, const double b1,
	const double default_x, const double default_y,
	const std::function<void(const double, const double)>& receiver
) {
	if (std::abs(b0) < TOL &&& std::abs(a0) < TOL && std::abs(c0) < TOL) {

		return;
	}
#if 0
	const double A = std::pow(a0 * b1, 2) + (a0 * c1) * (b0 * c1);
	const double B = 2 * (a0 * e1) * (a0 * b1) - (a0 * c1) * (c0 * e1) + (a0 * c1) * (b0 * d1) + (a0 * d1) * (b0 * c1);
	const double C = std::pow(a0 * e1, 2) + 2 * (a0 * b1) * (a0 * f1) - (a0 * c1) * (c0 * f1)
		- (a0 * c1) * (d0 * e1) - (a0 * d1) * (c0 * e1) + (a0 * d1) * (b0 * d1);
	const double D = 2 * (a0 * e1) * (a0 * f1) - (a0 * c1) * (d0 * f1) - (a0 * d1) * (c0 * f1) - (a0 * d1) * (d0 * e1);
	const double E = std::pow(a0 * f1, 2) - (a0 * d1) * (d0 * f1);
#endif
	const double q0 = -a0*a1*d0*d1*f0 + std::pow(a0, 2)*std::pow(d1, 2)*f0 + a0*std::pow(a1, 2)*std::pow(f0, 2) + a0*a1*std::pow(d0, 2)*f1 - std::pow(a0, 2)*d0*d1*f1 - 2*std::pow(a0, 2)*a1*f0*f1 + std::pow(a0, 3)*std::pow(f1, 2);
	const double q1 = -a0*a1*d0*d1*e0 + std::pow(a0, 2)*std::pow(d1, 2)*e0 + a0*a1*std::pow(d0, 2)*e1 - std::pow(a0, 2)*d0*d1*e1 - a0*a1*c1*d0*f0 - a0*a1*c0*d1*f0 + 2*std::pow(a0, 2)*c1*d1*f0 + 2*a0*std::pow(a1, 2)*e0*f0 - 2*std::pow(a0, 2)*a1*e1*f0 + 2*a0*a1*c0*d0*f1 - std::pow(a0, 2)*c1*d0*f1 - std::pow(a0, 2)*c0*d1*f1 - 2*std::pow(a0, 2)*a1*e0*f1 + 2*std::pow(a0, 3)*e1*f1;
	const double q2 = a0*a1*b1*std::pow(d0, 2) - a0*a1*b0*d0*d1 - std::pow(a0, 2)*b1*d0*d1 + std::pow(a0, 2)*b0*std::pow(d1, 2) - a0*a1*c1*d0*e0 - a0*a1*c0*d1*e0 + 2*std::pow(a0, 2)*c1*d1*e0 + a0*std::pow(a1, 2)*std::pow(e0, 2) + 2*a0*a1*c0*d0*e1 - std::pow(a0, 2)*c1*d0*e1 - std::pow(a0, 2)*c0*d1*e1 - 2*std::pow(a0, 2)*a1*e0*e1 + std::pow(a0, 3)*std::pow(e1, 2) + 2*a0*std::pow(a1, 2)*b0*f0 - 2*std::pow(a0, 2)*a1*b1*f0 - a0*a1*c0*c1*f0 + std::pow(a0, 2)*std::pow(c1, 2)*f0 - 2*std::pow(a0, 2)*a1*b0*f1 + 2*std::pow(a0, 3)*b1*f1 + a0*a1*std::pow(c0, 2)*f1 - std::pow(a0, 2)*c0*c1*f1;
	const double q3 = 2*a0*a1*b1*c0*d0 - a0*a1*b0*c1*d0 - std::pow(a0, 2)*b1*c1*d0 - a0*a1*b0*c0*d1 - std::pow(a0, 2)*b1*c0*d1 + 2*std::pow(a0, 2)*b0*c1*d1 + 2*a0*std::pow(a1, 2)*b0*e0 - 2*std::pow(a0, 2)*a1*b1*e0 - a0*a1*c0*c1*e0 + std::pow(a0, 2)*std::pow(c1, 2)*e0 - 2*std::pow(a0, 2)*a1*b0*e1 + 2*std::pow(a0, 3)*b1*e1 + a0*a1*std::pow(c0, 2)*e1 - std::pow(a0, 2)*c0*c1*e1;
	const double q4 = a0*std::pow(a1, 2)*std::pow(b0, 2) - 2*std::pow(a0, 2)*a1*b0*b1 + std::pow(a0, 3)*std::pow(b1, 2) + a0*a1*b1*std::pow(c0, 2) - a0*a1*b0*c0*c1 - std::pow(a0, 2)*b1*c0*c1 + std::pow(a0, 2)*b0*std::pow(c1, 2);

	std::cout << "q0=" << q0 << std::endl;
	std::cout << "q1=" << q1 << std::endl;
	std::cout << "q2=" << q2 << std::endl;
	std::cout << "q3=" << q3 << std::endl;
	std::cout << "q4=" << q4 << std::endl;

	if (std::abs(q0) < TOL && std::abs(q1) < TOL && std::abs(q2) << TOL && std::abs(q3) < TOL && std::abs(q4) < 0) {
		// TODO:
	}

	solve_4(
//		E, D, C, B, A,
		q0, q1, q2, q3, q4,
		default_y,
		[q0, q1, q2, q3, q4, a0, b0, c0, d0, e0, f0, a1, b1, c1, d1, e1, f1, default_x, &receiver](const double y) {
#if 0
			const double c0 = f0 + d0 * y + b0 * y * y;
			const double c1 = d0 + c0 * y;
			const double c2 = a0;
#endif
			std::cout << "\t\tquart 1: " << q0 + q1 * y + q2 * y * y + q3 * y * y * y + q4 * y * y * y * y << std::endl;

			const double l0 = a1*b0*std::pow(y, 2) - a0*b1*std::pow(y, 2) + a1*e0*y - a0*e1*y + a1*f0 - a0*f1;
			const double l1 = a1*c0*y - a0*c1*y + a1*d0 - a0*d1;

			const double pv0 = c0 * y + d0;
			const double qv0 = f0 + e0 * y + b0 * y * y;

			const double pv1 = c1 * y + d1;
			const double qv1 = f1 + e1 * y + b1 * y * y;

			std::cout << "l1=" << l1 << std::endl;
			std::cout << "l0=" << l0 << std::endl;
			std::cout << "p=" << pv0 << std::endl;
			std::cout << "q=" << qv0 << std::endl;
			std::cout << "a0=" << a0 << std::endl;

			const auto rec = [pv0, qv0, pv1, qv1, q0, q1, q2, q3, q4, a0, b0, c0, d0, e0, f0, a1, b1, c1, d1, e1, f1, l0, l1, y, &receiver](const double x){
				std::cout << "\t\tlinear: " << l0 + l1 * x << std::endl;
				std::cout << "\t\tquart 2: " << a0 * x * x + pv0 * x + qv0 << std::endl;
				std::cout << "\t\tquart 3: " << a1 * x * x + pv1 * x + qv1 << std::endl;
				std::cout << "\t\tinter linear: " << (
					a1*b0*y * y - a0*b1*y * y + a1*e0*y - a0*e1*y + a1*f0 - a0*f1 + (a1*c0*y - a0*c1*y + a1*d0 - a0*d1)*x
				) << std::endl;
				receiver(x, y);
			};

			if (std::abs(l1) > TOL) {
				rec(-l0 / l1);
				return;
			}
			solve_20_20(
				qv0, pv0, a0,
				qv1, pv1, a1,
				default_x,
				rec
			);
		}
	);
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

	if (0)
	{
		const double c0 = -0.034133333333333335;
		const double c1 = 0.17066666666666666;
		const double c2 = 0;
		const double c3 = -0.53333333333333333;
		// -0.034133333333333335 + 0.17066666666666666 * x + -0.53333333333333333 * x^2 + x^4

		solve_4n(
			c0, c1, c2, c3,
			[c0, c1, c2, c3](double x) {
				std::cout << "Quartic solution: " << x << std::endl;
				std::cout << "p(x) = " << (c0 + c1 * x + c2 * x * x + c3 * x * x * x + x * x * x * x) << std::endl;
			});
	}

	{
		const double a00 = -1;
		const double a10 = 0;
		const double a01 = 0;
		const double a20 = 1;
		const double a11 = 0;
		const double a02 = 1;

		const double b00 = -1;
		const double b10 = 0;
		const double b01 = 0;
		const double b20 = 1;
		const double b11 = 0;
		const double b02 = 0;
		solve_22(
			a00, a10, a01, a20, a11, a02,
			b00, b10, b01, b20, b11, b02,
			std::nan(""), std::nan(""),
			[
				a00, a10, a01, a20, a11, a02,
				b00, b10, b01, b20, b11, b02](const double x, const double y
			) {
				std::cout << "=========================================================" << std::endl;
				std::cout << "systems solution:" << std::endl;
				std::cout << x << ", " << y << std::endl;

				std::cout << "eq1: " << a00 + a10 * x + a01 * y + a20 * x * x + a11 * x * y + a02 * y * y << std::endl;
				std::cout << "eq2: " << b00 + b10 * x + b01 * y + b20 * x * x + b11 * x * y + b02 * y * y << std::endl;
				std::cout << "=========================================================" << std::endl;
			}
		);

		{
			const double x = 1;
			const double y = 0;
			std::cout << "=========================================================" << std::endl;
			std::cout << "some solution" << std::endl;
			std::cout << "eq1: " << a00 + a10 * x + a01 * y + a20 * x * x + a11 * x * y + a02 * y * y << std::endl;
			std::cout << "eq2: " << b00 + b10 * x + b01 * y + b20 * x * x + b11 * x * y + b02 * y * y << std::endl;
			std::cout << "=========================================================" << std::endl;
		}
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