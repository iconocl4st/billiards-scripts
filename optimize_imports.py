#!/bin/python3

import re
import os


PROJECTS = [
	'billiards-attempts-api',
	'billiards-common',
	'billiards-graphics-api',
	'billiards-layouts-api',
	#'billiards-scripts',
	'billiards-client',
	'billiards-config-api',
	'billiards-image-processing-api',
	'billiards-projection-api',
	'billiards-shots-api',
]

PROJECT_INCLUDES = re.compile('#include\s*"([^"]*)"')
SYSTEM_PATTERN = re.compile('#include\s*<([^>]*)>')


class Headers:
	def __init__(self):
		self.system_includes = set()
		self.project_includes = set()
	
	def scan(self, line):
		m = PROJECT_INCLUDES.search(line)
		if (m):
			self.project_includes.add(m.group(1))
			
		m = SYSTEM_PATTERN.search(line)
		if (m):
			self.system_includes.add(m.group(1))


def get_repos_dir():
	if 'REPOS' in os.environ and os.path.exists(os.environ['REPOS']):
		repos_dir = os.environ['REPOS']
	else:
		repos_dir = os.path.join(__file__, '..', '..')
	return os.path.realpath(repos_dir)


def parse_header(path):
	headers = Headers()
	with open(path, 'r') as infile:
		for line in infile:
			headers.scan(line)
	return headers


def include_file(filename):
	if filename[-len('.h'):] == '.h':
		return True
	if filename[-len('.cpp'):] == '.cpp':
		return True
	if filename[-len('.c'):] == '.c':
		return True
	return False


def collect_headers(project, headers):
	repos = get_repos_dir()
	project_root = os.path.join(repos, project, 'src')
	
	for root, dirs, files in os.walk(project_root):
		for filename in files:
			if not include_file(filename):
				continue
			
			relative_path = root[len(project_root):]
			abs_path = os.path.join(root, filename)
			headers[(project_root, os.path.join(relative_path, filename))] = parse_header(abs_path)


def locate_header(headers, include):
	for key, value in headers.items():
		if key[1] == include:
			return value
	return None


def print_headers(headers, current, depth):
	indents = ''.join(['\t'] * depth)
	
	child_includes = set()
	current_includes = set()
	
	print(indents + 'system includes:')
	for system_include in current.system_includes:
		print(indents + '\t' + system_include)
		current_includes.add('<' + system_include + '>')
		
	print(indents + 'project includes:')
	for project_include in current.project_includes:
		print(indents + '\t' + project_include)
		current_includes.add('"' + project_include + '"')
		header = locate_header(headers, project_include)
		if header is not None:
			included = print_headers(headers, header, depth + 2)
			current_includes.update(included)
	
	print(indents + 'redundant:')
	for include in current_includes:
		if include in child_includes:
			print(indents + '\t' + include)
	
	r = set()
	r.update(child_includes)
	r.update(current_includes)
	return r


def main():
	headers = {}
	for project in PROJECTS:
		collect_headers(project, headers)
	for key, value in headers.items():
		if key[1][-len('.cpp'):] == '.cpp':
			print(key)
			print_headers(headers, value, 1)
	


if __name__ == '__main__':
	main()


