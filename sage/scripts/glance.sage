

'''

load("glance.sage"); determine_glance_location()


'''

# Rolling glance direction.
# Using the model in https://billiards.colostate.edu/technical_proofs/new/TP_A-4.pdf

replacements = [
	PowerReplacement('dx', 2),
	PowerReplacement('dy', 2),
	PowerReplacement('ax', 4),
	PowerReplacement('ay', 4),
	PowerReplacement('r', 6),
]

def replaces(s):
	for replacement in replacements:
		s = replacement.replace(s)
	return s

def determine_glance_location():
	for replacement in replacements:
		replacement.generate_assignments()
	
	# The object ball is centered at the origin
	P.<x, y> = PolynomialRing(QQ, 2)
	# x, y = var('x y') # Location of the glance
	dx, dy = var('dx dy') # Location of the desired cue destination
	ax, ay = var('ax ay') # Location of the current cue position
	r = var('r') # The sum of radaii of the cue and object balls
	alpha = 2/7
	#if True:
	#	ax = -20; ay = 20
	#	dx = 17.759757877423255; dy = -4.7018134271419401
	#	r = 2.26

	# The tangent direction must be one of the following. Note, it will have magnitude r.
	for tdir_x, tdir_y in [(y, -x), (-y, x)]:
		# We need the location of the glance to be at radius r from the object ball
		poly_r = x^2 + y^2 - r^2
		
		# There is a point on the tangent line
		tx = x + tdir_x
		ty = y + tdir_y

		# Let s represent the distance travelled along the aim line such it is
		# orthogonal to the tangent line:
		s = var('s')
		aim_x = x + s * (x - ax); aim_y = y + s * (y - ay);
		poly_s = (aim_x - x) * (aim_x - tx) + (aim_y - y) * (aim_y - ty)
		poly_s = poly_s / s

		# The rolling glance means that the direction from the glance point satisfies
		t = var('t')
		eq_dx = t * (dx - x) == alpha * (tx - x) + (1 - alpha) * (aim_x - x)
		eq_dy = t * (dy - y) == alpha * (ty - y) + (1 - alpha) * (aim_y - y)
		eq_d = eq_dx * eq_dy.left() / t - eq_dy * eq_dx.left() / t
		poly_d = (eq_d - eq_d.right()).left()
		
		s_sol = solve(poly_s == 0, s)[0].right()
		poly_d = poly_d.substitute(s=s_sol)
		poly_d = (poly_d * s_sol.denominator()).simplify_rational()

		res = poly_d.resultant(poly_r, y)
		#res = (7.105427357601002e-15)*x^6 + res
		#res_roots = res.roots()
		#for r in res_roots:
		#	numerical_approx(r[0])
			
		
		print("\t\t{")
		for coef, ord in res.coefficients(x):
			print_double_assignment('a' + str(ord), coef, replaces=replaces)
		print("\t\t\tstd::vector<double> coefficients{" + ", ".join('a' + str(ord) for _, ord in res.coefficients(x)) + "};")
		print_double_assignment('tx', tx)
		print_double_assignment('ty', ty)
		print_double_assignment('s_numerator', s_sol.numerator(), replaces=replaces)
		print_double_assignment('s_denominator', s_sol.denominator(), replaces=replaces)
		print('\t\t\tif (std::abs(s_denominator) < TOLERANCE) {\n\t\t\t\treturn;\n\t\t\t}')
		print('\t\t\tconst double s = s_numerator / s_denominator;')
		print_double_assignment('aim_x', aim_x)
		print_double_assignment('aim_y', aim_y)
		print('\t\t\tchecker(RollingGlanceCalculation{x, y, tx, ty, aim_x, aim_y});')
		print("\t\t}")
		
		
		#for coef, ord in poly_s.coefficients(s):
		#	print_double_assignment('b' + str(ord), coef, replaces=replaces)
		#print("\t\t\tstd::vector<double> s_coefficients{" + ", ".join('b' + str(ord) for _, ord in poly_s.coefficients(s)) + "};")
		
		
		# print(solve([poly_r == 0, poly_d == 0], x, y))
		# We need to solve for x and y
		# First, we solve for s. One solution is s = 0, so we divide it out.

		
		#poly_d = (poly_d * s_sol.denominator()).simplify_rational()
		#poly_d = poly_d.expand()
		
		# print(res.factor())
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
