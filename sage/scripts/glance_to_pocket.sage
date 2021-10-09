
# export REPOS=/mnt/1f0ab4b3-c472-49e1-92d8-c0b5664f7fdb/ProjectsForFun/Pool/repos
# docker run -v $(realpath $REPOS/billiards-scripts/sage/scripts):/home/sage/scripts --workdir=/home/sage/scripts -it sagemath/sagemath:latest

from .common import convert_powers, get_coefficients

def determine_dest_signs():
	# The desired cue location is along a line that is orthogonal to the pocket at some point dx, dy
	# eq_d = (dx - x) * (dx - px) + (dy - y) * (dy - py) == 0
	dx, dy = var('dx dy')
	dx_eqns = [
		[dx - x == dy - py, dy - y == -(dx - px)],
		[dx - x == -(dy - py), dy - y == dx - px],
	]
	# px * x + py * y >= 0
	# px^2 + py^2 >= 0
	for dxdy_eqns in dx_eqns:
		dxdy_sols = solve(dxdy_eqns, dx, dy)[0]
		print(dxdy_sols)
		dx_sol = dxdy_sols[0].right()
		dy_sol = dxdy_sols[1].right()
		print((dx_sol * px + dy_sol * py).factor() >= 0)

	#if px * y - py * x > 0
	#   tx = y, ty = -x
	#
	# else
	#   tx = -y, tx = x
#determine_dest_signs()


def determine_dest_signs_as_direction():
	# The desired cue location is along a line that is orthogonal to the pocket at some point dx, dy
	# eq_d = dx * (x + dx - px) + dy * (y + dy - py) == 0
	x, y, dx, dy, px, py = var('x y dx dy px py')
	dx_eqns = [
		[dx == y + dy - py, dy == -(x + dx - px)],
		[dx - x == -(y + dy - py), dy - y ==x + dx - px],
	]
	# px * x + py * y >= 0
	# px^2 + py^2 >= 0
	for dxdy_eqns in dx_eqns:
		dxdy_sols = solve(dxdy_eqns, dx, dy)[0]
		print(dxdy_sols)
		dx_sol = dxdy_sols[0].right()
		dy_sol = dxdy_sols[1].right()
		print(((dx_sol - x)* px + (dy_sol - y) * py).factor() >= 0)

	#if px * y - py * x > 0
	#   tx = y, ty = -x
	#   1/2*px - 1/2*py - 1/2*x + 1/2*y, dy == 1/2*px + 1/2*py - 1/2*x - 1/2*y
	# else
	#   tx = -y, tx = x
# determine_dest_signs_as_direction()


def replace_vars(e):
	return (e
		.replace('std::pow(r1, 2)', 'r12')
		.replace('std::pow(r1, 3)', 'r13')
		.replace('std::pow(r1, 4)', 'r14')
		.replace('std::pow(r2, 2)', 'r22')
		.replace('std::pow(r2, 3)', 'r23')
		.replace('std::pow(r2, 4)', 'r24')
		.replace('std::pow(px, 2)', 'px2')
		.replace('std::pow(px, 3)', 'px3')
		.replace('std::pow(px, 4)', 'px4')
		.replace('std::pow(py, 2)', 'py2')
		.replace('std::pow(py, 3)', 'py3')
		.replace('std::pow(py, 4)', 'py4')
		#.replace('std::pow(ax, 2)', 'ax2')
		#.replace('std::pow(ay, 2)', 'ay2')
		#.replace('px2 + py2', 'px2ppy2')
		#.replace('ax2 + ay2', 'ax2pay2')
		#.replace('sqrt(px2ppy2 - r12)', 'sradi')
		#.replace('(double) (px2ppy2)', 'px2ppy2')
		#.replace('(double) (ax2pay2)', 'ax2pay2')
		#.replace('(double) py', 'py')
		#.replace('(double) px', 'px')
	)

# load("glance_to_pocket.sage"); rolling_glance_to_pocket()

def rolling_glance_to_pocket():
	# Rolling glance direction.
	# Using the model in https://billiards.colostate.edu/technical_proofs/new/TP_A-4.pdf

	P.<x, y> = PolynomialRing(QQ, 2)
	# The object ball is centered at the origin
	ax, ay = var('ax ay') # Location of the current cue position
	px, py = var('px py') # Location of the pocket
	r1 = var('r1') # The sum of radius of the cue ball
	r2 = var('r2') # The sum of radius of the object ball
	alpha = 2/7


	# The tangent direction must be orthogonal to (x, y)
	for is_valid, tx, ty, dx, dy in [(
			'px * y - py * x > 0',
			y, -x,
			1/2*px - 1/2*py - 1/2*x + 1/2*y, 1/2*px + 1/2*py - 1/2*x - 1/2*y
		), (
			'px * y - py * x < 0',
			-y, x,
			1/2*px + 1/2*py - y, -1/2*px + 1/2*py + x
	)]:
		# When pocket is on the same sid
		print('\t\t// if ' + is_valid)
		print('\t\t// tangent direction (tx, ty) = ' + str([tx, ty]))
		print('\t\t// destination direction (dx, dy) = (' + str(dx) + ',' + str(dy) + ')')

		# We need the location of the glance to be at radius r1 + r2 from the object ball
		eq_r1 = x^2 + y^2 == (r1 + r2)^2
		poly_r1 = (eq_r1 - eq_r1.right()).left()

		# Let s1 be a distance travelled along the aim line such that it is orthogonal to the tangent line:
		s1 = var('s1')
		eq_orth = ax * (s1 * ax - tx) + ay * (s1 * ay - ty) == 0
		s1 = solve(eq_orth, s1)[0].right()
		print('\t\t// computed orthogonal requirement: s1=' + str(s1))

		# The rolling glance means that the direction from the glance point satisfies
		s3 = var('s3')
		eq_dx = s3 * dx == (1 - alpha) * tx + alpha * s1 * ax
		eq_dy = s3 * dy == (1 - alpha) * ty + alpha * s1 * ay

		eq_s = dy * eq_dx - dx * eq_dy
		poly_s = (eq_s - eq_s.left()).right()
		poly_s = (poly_s * s1.denominator()).simplify_rational()

		# s_sol = solve(eq_dx, s)[0].right()
		# eq_d = (eq_dy.substitute(s=s_sol) * s_sol.denominator()).simplify_rational().expand()
		# poly_d = (eq_d - eq_d.right()).left().factor()
		#print('\t\t// computed rolling requirement: poly_d=' + str(poly_d))

		# The desired location dx, dy we must be at distance r1 from the pocket
		#s4 = var('s4')
		eq_r2 = (x + s4 * dx - px)^2 + (y + s4 * dy - py)^2 == r1^2
		poly_r2 = (eq_r2 - eq_r2.right()).left()

		# Variables: x, y
		# Equations: poly_r1 == 0, poly_r2 == 0, poly_d == 0,
		print(poly_r1.expand())
		print(poly_r2.expand())
		print(poly_s.expand().factor())


		print('\t\t{')
		for coef, name in get_coefficients(poly_r1.expand(), x, y):
			print('\t\t\tconst double a' + str(name) + ' = ' + replace_vars(convert_powers(coef)) + ';')
		for coef, name in get_coefficients(poly_s.expand(), x, y):
			print('\t\t\tconst double b' + str(name) + ' = ' + replace_vars(convert_powers(coef)) + ';')
		print('\t\t\tconst auto get_glance = [&](const double x, const double y) {');
		print('\t\t\t\treturn KissToPocketSolution{')
		print('\t\t\t\t\tx, y,')
		print('\t\t\t\t\t' + convert_powers(str(dx)) + ', ' + convert_powers(str(dy)) + ',')
		print('\t\t\t\t\t' + convert_powers(str(tx)) + ', ' + convert_powers(str(ty)) + ',')
		print('\t\t\t\t\tax - x, ay - y')
		print('\t\t\t\t};')
		print('\t\t\t};')
		print('\t\t}')
		if True:
			continue


		# r1 = poly_r1.resultant(poly_s, y).expand().factor()
		#print("poly")
		#print(r1.factor())
		#print(r1.operands()[0].operands()[0].factor())
		# solve([poly_r1, poly_r2, poly_s], x, y, s3)
		if True:
			continue

		r1 = poly_r1.resultant(poly_r2, y)
		print('\t\t// first resultant polynomial: ' + str(r1))
		r2 = r1.resultant(poly_d, y)
		print('\t\t// second resultant polynomial: ' + str(r2))
		poly = r2.factor()
		print('\t\t// final resultant polynomial: ' + str(poly))

		factor = poly.operands()[0].operands()[0]
		if not (factor^2 / 16 - r2).is_zero():
			raise Exception('assumption error')

		print('\t\t{')
		for coef, ord in factor.expand().coefficients(x):
			print('\t\t\tconst double c' + str(ord) + ' = ' + replace_vars(convert_powers(coef)) + ';')
		print('\t\t}')
