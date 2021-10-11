
import re


# export REPOS=/mnt/1f0ab4b3-c472-49e1-92d8-c0b5664f7fdb/ProjectsForFun/Pool/repos
# docker run -v $(realpath $REPOS/billiards-scripts/sage/scripts):/home/sage/scripts --workdir=/home/sage/scripts -it sagemath/sagemath:latest



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


class PowerReplacement:
	def __init__(self, term, power):
		self.term = term
		self.power = power
	
	def get_var_name(self, p):
		if p == 1:
			return self.term
		return self.term + "_" + str(p)
	
	def replace(self, s):
		for p in range(2, self.power + 1):
			s = s.replace('std::pow(' + self.term + ', ' + str(p) + ')', self.get_var_name(p))
		return s
	
	def generate_assignments(self, indent=2):
		for p in range(2, self.power + 1):
			print_double_assignment(self.get_var_name(p), self.get_var_name(p - 1) + ' * ' + self.term);

# p = PowerReplacement('ax', 4)
# p.generate_assignments()
# print(p.replace('std::pow(ax, 3)'))
