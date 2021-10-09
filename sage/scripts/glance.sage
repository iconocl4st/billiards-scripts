

'''



load("glance.sage"); determine_glance_location()


'''

# from common import *

# Rolling glance direction.
# Using the model in https://billiards.colostate.edu/technical_proofs/new/TP_A-4.pdf

def replaces(s):
	return (s
			.replace('std::pow(dx, 2)', 'dx2')
			.replace('std::pow(dy, 2)', 'dy2')
			.replace('std::pow(ax, 2)', 'ax2')
			.replace('std::pow(ax, 3)', 'ax3')
			.replace('std::pow(ax, 4)', 'ax4')
			.replace('std::pow(ay, 2)', 'ay2')
			.replace('std::pow(ay, 3)', 'ay3')
			.replace('std::pow(ay, 4)', 'ay4')
			.replace('std::pow(r, 2)', 'r2')
	)

def determine_glance_location():
	P.<x, y> = PolynomialRing(QQ, 2)
	# The object ball is centered at the origin
	dx, dy = var('dx dy') # Location of the desired cue destination
	ax, ay = var('ax ay') # Location of the current cue position
	r = var('r') # The sum of radaii of the cue and object balls
	alpha = 2/7

	# The tangent direction must be one of the following. Note, it will have magnitude r.
	for tx, ty in [(y, -x), (-y, x)]:
		# We need the location of the glance to be at radius r from the object ball
		eq_r = x^2 + y^2 == r^2
		poly_r = (eq_r - eq_r.right()).left()

		# Let s represent the distance travelled along the aim line such it is
		# orthogonal to the tangent line:
		s = var('s')
		aim_x = x - ax; aim_y = y - ay;
		eq_s = aim_x * (s * aim_x - tx) + ay * (s * aim_y - ty) == 0
		s_sol = solve(eq_s, s)[0].right()

		# The rolling glance means that the direction from the glance point satisfies
		t = var('t')
		eq_dx = t * dx == alpha * tx + (1 - alpha) * s_sol * aim_x
		eq_dy = t * dy == alpha * ty + (1 - alpha) * s_sol * aim_y
		eq_d = dy * eq_dx - dx * eq_dy
		poly_d = (eq_d - eq_d.right()).left()
		poly_d = (poly_d * s_sol.denominator()).simplify_rational()
		poly_d = poly_d.expand()

		# We need to solve for x and y
		# First, we solve for s. One solution is s = 0, so we divide it out.

		res = poly_r.resultant(poly_d, y)
		print(res.factor())
		#print(res.expand().coefficients(x))
		# print(poly_d)
		#
		# print("\t\t{")
		# for coef, idx in get_coefficients(poly_r, x, y):
		# 	print_double_assignment('a' + idx, coef, replaces=None)
		# for coef, idx in get_coefficients(poly_d, x, y):
		# 	print_double_assignment('b' + idx, coef, replaces=None)
		#
		# print_double_assignment('tx', tx)
		# print_double_assignment('ty', ty)
		# print('\t\t\tchecker(RollingGlanceCalculation{x, y, tx, ty, x - ax, y - ay});')
		# print("\t\t}")
