
		#expr = eq_dx / eq_dy; expr
		#print(expr)
		if True:
			return
		expr = expr * ((alpha*ay*s - (alpha - 1)*ty)); expr
		expr = expr - (alpha*ax*s - (alpha - 1)*tx); expr
		expr = expr.left(); expr
		expr = expr * dy; expr
		expr = expr.expand(); expr
		expr = expr.substitute(s=s_sol); expr
		expr = expr * (ax^2 + ay^2); expr
		expr = expr.simplify_rational(); expr
		expr = expr.expand(); expr
		expr = expr.substitute(alpha=2/7); expr

		print('\t\t// tangent assignment:')
		print('\t\tconst double tx  = ' + str(tx) + ';')
		print('\t\tconst double ty  = ' + str(ty) + ';')

		y_sol = solve(expr.substitute(**tan_sol), y)[0].right()
		expr2 = eq_r.substitute(y=y_sol)

		expr2 = expr2 - r^2
		expr2 = expr2 * (2*ax*ay*dx - (7*ax^2 + 5*ay^2)*dy)^2
		expr2 = expr2.simplify_rational()
		expr2 = expr2.left()
		expr2 = expr2.expand()
		expr2 = expr2.collect(x)
		for c, ord in expr2.coefficients(x):
			print('\t\tconst double c' + str(ord) + ' = ' + convert_powers(c.simplify())
				.replace('std::pow(ax, 2)', 'ax2')
				.replace('std::pow(ax, 3)', 'ax3')
				.replace('std::pow(ax, 4)', 'ax4')
				.replace('std::pow(ay, 2)', 'ay2')
				.replace('std::pow(ay, 3)', 'ay3')
				.replace('std::pow(ay, 4)', 'ay4')
				.replace('std::pow(dx, 2)', 'dx2')
				.replace('std::pow(dy, 2)', 'dy2')
				.replace('std::pow(r, 2)', 'r2')
				+ ';')

				
				
				
				
				
				
				
				
				
				
				
