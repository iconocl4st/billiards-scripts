
def determine_dest_signs():
	# The desired cue location is along a line that is orthogonal to the pocket at some point dx, dy
	# eq_d = (dx - x) * (dx - px) + (dy - y) * (dy - py) == 0
	x, y = var('x y')
	px, py = var('px py')
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
# determine_dest_signs()


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
#determine_dest_signs_as_direction()

replacements = [
	PowerReplacement('r1', 8),
	PowerReplacement('r2', 8),
	PowerReplacement('ax', 4),
	PowerReplacement('ay', 4),
	PowerReplacement('px', 2),
	PowerReplacement('py', 2),
]


def replaces(s):
	for replacement in replacements:
		s = replacement.replace(s)
	return s
	
def use_solution_2(poly, var, solution):
	poly = poly.expand()
	deg = poly.degree(var)
	num = solution.numerator()
	den = solution.denominator()
	ret = 0
	for coefficent, order in poly.coefficients(var):
		ret = ret + coefficent * num^order * den^(deg - order)
	ret = ret.expand()
	return ret

def use_solution(poly, var, solution):
	poly = poly.substitute(**{var: solution})
	poly = poly * solution.denominator()
	poly = poly.simplify_rational()
	poly = poly.expand()
	return poly

	
class Linearization:
	def __init__(self):
		self.mapping = {}
		self.rmap = []
		self.index = 0
	
	def get_mapping(self, order):
		if order not in self.mapping:
			self.mapping[order] = self.index
			self.rmap.append(order)
			self.index += 1
		return self.mapping[order]
	
	def print(self, outf):
		outf.write('\t\tstd::vector<std::string> mapping {' + 
				', '.join('"' + str(pow)+ '"' for pow in self.rmap)
			+ '};\n')
	
def print_assign_to_matrix(poly, row, vars, lin, outf):
	outf.write('\t\t//' + str(poly) + '\n')
	for term in poly.expand().operands():
		powers = tuple(term.degree(var) for var in vars)
		expr = term / prod(var^p for var, p in zip(vars, powers))
		outf.write('\t\t// ' + str(powers) + ': ' + str(term) + '\n')
		outf.write('\t\tMat(' + str(row) + ', ' + str(lin.get_mapping(powers)) + ') += ' + replaces(convert_powers(expr)) + ';\n')

'''

load("glance_to_pocket.sage");

with open('./output.c', 'w') as outf:
	rolling_glance_to_pocket_2(outf)

'''


def rolling_glance_to_pocket(outf, idx):
	# Rolling glance direction.
	# Using the model in https://billiards.colostate.edu/technical_proofs/new/TP_A-4.pdf
	
	for replacement in replacements:
		replacement.generate_assignments(outf=outf)
	
	# P.<x, y, u> = PolynomialRing(QQ, 3)
	x, y = var('x y') # Location of glance
	# The object ball is centered at the origin
	ax, ay = var('ax ay') # Location of the current cue position
	px, py = var('px py') # Location of the pocket
	r1 = var('r1') # The radius of the cue ball
	r2 = var('r2') # The radius of the object ball
	alpha = 2/7


	# The tangent direction must be orthogonal to (x, y)
	is_valid, tdir_x, tdir_y, ddir_x, ddir_y = [(
			'px * y - py * x > 0',
			y, -x,
			1/2*px - 1/2*py - 1/2*x + 1/2*y, 1/2*px + 1/2*py - 1/2*x - 1/2*y
		), (
			'px * y - py * x < 0',
			-y, x,
			1/2*px + 1/2*py - y, -1/2*px + 1/2*py + x
	)][idx]
	# We need the location of the glance to be at radius r from the object ball
	poly_r1 = x^2 + y^2 - (r1 + r2)^2
	poly_r1 = poly_r1.expand()
	
	# There is a point on the tangent line
	tx = x + tdir_x
	ty = y + tdir_y
	
	u = var('u')
	dx = x + u * ddir_x
	dy = y + u * ddir_y

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
	poly_d = (14 * poly_d / u).simplify_rational()
		
	# The desired location dx, dy we must be at distance r1 from the pocket
		
	poly_r2 = (dx - px)^2 + (dy - py)^2 - r1^2
	poly_r2 = poly_r2.expand()
		
	# The desired location must be orthogonal to the pocket
	# poly_n = (dx - x) * (dx - px) + (dy - y) * (dy - py)
	# poly_n = poly_n / u
		
	s_sol = solve(poly_s == 0, s)[0].right()
	poly_d = use_solution(poly_d, 's', s_sol)
	
	# Variables: x, y, u
	# Equations: poly_r1 == 0, poly_r2 == 0, poly_d == 0
	outf.write('====\n')
	
	print('r1')
	print(poly_r1)
	print('r2')
	print(poly_r2)
	print('d')
	print(poly_d)
	
	poly_d = poly_d.maxima_methods().divide(poly_r1)[1].expand()
	#poly_d = poly_d.maxima_methods().divide(poly_r2)[1]
	#poly_r2 = poly_r2.maxima_methods().divide(poly_r1)[1]
	#poly_r2 = poly_r2.maxima_methods().divide(poly_d)[1]
	
	#print('division:')
	#for idx, r in enumerate(r1):
	#	print(idx, r.expand())
	#print((r1[0] * poly_r1 + r1[1] - poly_d).is_zero())
	
	#print('r1')
	#print(poly_r1)
	#print('r2')
	#print(poly_r2)
	print('d')
	print(poly_d.resultant(poly_r1, y))
	print(solve(poly_d == 0, y))
	
	if True:
		return
		
		
	r = poly_r1.resultant(poly_d, y)
	print(r.expand())
		
		
	outf.write('\t\t{\n')
	for coef, ord in r.expand().coefficients(x):
		print_double_assignment('a' + str(ord), coef, replaces=replaces, outf=outf)
	outf.write('\t\t}\n')
		
		
	outf.write('\t\t{\n')
	for coef, ord in poly_r2.expand().coefficients(u):
		print_double_assignment('b' + str(ord), coef, replaces=replaces, outf=outf)
	outf.write('\t\t}\n')
	
	if True:
		return
		
		
		
		
		
		
		
		
		
		
		
		
	vars = [x, y, u]
	lin = Linearization()
	outf.write('\t\t{\n')
	print_assign_to_matrix(poly_r1, 0, vars, lin, outf)
	print_assign_to_matrix(poly_r2, 1, vars, lin, outf)
	print_assign_to_matrix(poly_d, 2, vars, lin, outf)
	outf.write('\t\t\tEigen::MatrixXd Mat{3, ' + str(lin.index) + '};\n')
	outf.write('\t\t\tfor (int i=0;i<3;i++) for (int j=0;j<' + str(lin.index) + ';j++) Mat(i, j) = 0;\n')
		
	outf.write('\t\t\tstd::cout << "Found matrix: " << std::endl;\n')
	outf.write('\t\t\tstd::cout << Mat << std::endl;\n')
	outf.write('\t\t}\n')
	lin.print(outf)
	# print('=================================')
		
		
		
		
		
		
		
		
	#poly = poly_d
	#poly = poly - poly.coefficient(y^4) * y^2 * poly_r1
	#poly = poly.expand()
	#poly = poly - poly.coefficient(y^3) * y * poly_r1
	#poly = poly.expand()
	#poly = poly - poly.coefficient(y^2) * poly_r1
	#poly = poly.expand()
	#print(replaces(str(poly)))
		
	#y_sol = solve(poly == 0, y)[0].right()
	#poly_d = use_solution_2(poly_d, y, y_sol)
		
	#print(replaces(str(poly_d)))
	if True:
		return
		
	outf.write('\t\t{\n')
	for coef, ord in poly_d.coefficients(x):
		#for addend in coef.operands():
		#	for factor in addend.operands():
		#		print(factor.operator())
		#		print(factor)
		# print(replaces(str(coef.factor())))
		print_double_assignment('a' + str(ord), coef, replaces=replaces, outf=outf)
	#print_double_assignment('
	outf.write('\t\t}\n')
	lin.print()
	
	if True:
		return
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	#print([c[0].operands()[0] for c in poly_d.coefficients(x)])
			
	#print(poly_r2.coefficients())
		
		
	#J = P.ideal(poly_r1, poly_r2, poly_d)
	#print(J.dimension())

	#r1 = poly_r2.resultant(poly_d, u)
	#print(r1)
		
	#r1 = poly_r1.resultant(poly_r2, y)
	#
	#print(r1)
	#print(r1)
	#print(r2)
	#print('==')
	

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

	# r1 = poly_r1.resultant(poly_s, y).expand().factor()
	#print("poly")
	#print(r1.factor())
	#print(r1.operands()[0].operands()[0].factor())
	# solve([poly_r1, poly_r2, poly_s], x, y, s3)
	if True:
		return

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
	
		
		

def rolling_glance_to_pocket_2(outf):
	rolling_glance_to_pocket(outf, 0)
	rolling_glance_to_pocket(outf, 1)
