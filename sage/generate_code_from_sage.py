# docker run -v `pwd`:/home/sage/output -it sagemath/sagemath:latest


import json
from sage.symbolic.operators import add_vararg, mul_vararg, is_SymbolicVariable

def get_expr_type(expr):
	op = expr.operator()
	if op == mul_vararg:
		return '*'
	elif op == add_vararg:
		return '+'
	elif '<built-in function eq>' in str(op):
		return '=='
	elif '<built-in function pow>' in str(op):
		if len(expr.operands()) != 2:
			import pdb; pdb.set_trace()
		return '**'
	elif op is None:
		import pdb; pdb.set_trace()
		raise Exception()
	else:
		import pdb; pdb.set_trace()
		raise Exception()
	# return "Unknown: " + str(op) + ", " + str(type(op))


def get_expr(expr):
	if is_SymbolicVariable(expr):
		return {
			'type': 'variable',
			'name': str(expr)
		}

	if expr.operator() is None:
		return {
			'type': 'value',
			'value': str(expr)
		}
	return {
		'type': 'operator',
		'operator-type': get_expr_type(expr),
		'operands': [
			get_expr(operand)
			for operand in expr.operands()
		]
	}

class State:
	def __init__(self):
		self.commands = []
		self.var_count = 0
		self.args = set()

	def print(self, last_var, num, out):
		out.write('[[nodiscard]] inline' + '\n')
		out.write('bool get_solution_' + str(num).zfill(3) + '(' + '\n')
		for arg in self.args:
			out.write('\t\tconst double ' + arg + ',' + '\n')
		out.write('\t\tdouble *out_param' + '\n')
		out.write(') {' + '\n')
		for cmd in self.commands:
			out.write(cmd + '\n')
		out.write('\t*out_param = ' + last_var + ';' + '\n')
		out.write('\treturn true;' + '\n')
		out.write('}' + '\n')

	def next_var(self):
		num = self.var_count
		self.var_count += 1
		return 'var_' + str(num).zfill(3)


def expr_to_c(j, state):
	if j['type'] == 'operator':
		if len(j['operands']) < 2:
			print('here')
		prev = None
		for operand in j['operands']:
			var_name = expr_to_c(operand, state)

			if prev is not None:
				next_var = state.next_var()
				decl = 'const double ' + next_var + ' = '

				guard = None
				e = 'todo'
				if j['operator-type'] == '*':
					e = prev + ' * ' + var_name
				if j['operator-type'] == '+':
					e = prev + ' + ' + var_name
				elif j['operator-type'] == '**':
					if var_name == '2':
						e = prev + ' * ' + prev
					elif var_name == '0.5':
						guard = 'if (' + prev + ' < 0) { return false; }'
						e = 'std::sqrt(' + prev + ')'
					elif var_name == '-1':
						guard = 'if (std::abs(' + prev + ') < TOLERANCE) { return false; }'
						e = '1.0 / ' + prev
					else:
						e = 'std::pow(' + prev + ', ' + var_name + ')'
				elif j['operator-type'] == '==':
					raise Exception()

				if guard is not None:
					state.commands.append('\t' + guard)
				state.commands.append('\t' + decl + e + ';')
				prev = next_var
			else:
				prev = var_name
		return prev
	elif j['type'] == 'value':
		if j['value'] == '1/2':
			return '0.5';
		return j['value'].replace('/', '/ (double) ')
	elif j['type'] == 'variable':
		state.args.add(j['name'])
		return j['name']


def to_c(j, num, out):
	state = State()
	last_var = expr_to_c(j, state)
	state.print(last_var, num, out)


def print_sols(sols):
	with open('/home/sage/output/gen.cpp', 'w') as out:
		for idx, sol in enumerate(sols):
			j = get_expr(sol.right())
			to_c(j, idx, out)
		# print(json.dumps(j, indent=2))

if sols is not None:
	print_sols(sols)

