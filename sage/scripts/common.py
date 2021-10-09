
import re


# TODO: should use std::sqrt and std::cbrt...
def convert_powers(expr):
	return re.sub(
		r'([a-zA-Z0-9]*)\^([0-9]+)', r'std::pow(\g<1>, \g<2>)', str(expr)
	).replace('/', '/ (double) ')

def get_coefficients(expr, x, y):
	c00 = 0; c10 = 0; c01 = 0; c20 = 0; c11 = 0; c02 = 0;
	for term in expr.expand().operands():
		if   term.degree(x) == 0 and term.degree(y) == 0:
			c00 = c00 + term
		elif term.degree(x) == 1 and term.degree(y) == 0:
			c10 = c10 + term.coefficient(x)
		elif term.degree(x) == 0 and term.degree(y) == 1:
			c01 = c01 + term.coefficient(y)
		elif term.degree(x) == 2 and term.degree(y) == 0:
			c20 = c20 + term.coefficient(x^2)
		elif term.degree(x) == 1 and term.degree(y) == 1:
			c11 = c11 + term.coefficient(x*y)
		elif term.degree(x) == 0 and term.degree(y) == 2:
			c02 = c02 + term.coefficient(y^2)
		else:
			print(term)
			raise Exception('higher orders not supported')
	ret = [c00, c10, c01, c20, c11, c02]
	if not (c00 + c10 * x + c01 * y + c20 * x^2 + c11 * x * y + c02 * y^2 - expr).is_zero():
		raise Exception('Failure to get coefficients')
	return zip(ret, ['00', '10', '01', '20', '11', '02'])


def print_double_assignment(name, expression, indent=3, replaces=None):
	s = ''.join(['\t'] * indent + ['const double ', name, ' = ', convert_powers(expression), ';'])
	if replaces is not None:
		s = replaces(s)
	print(s)
