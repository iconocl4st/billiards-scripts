
from git import Repo
import sys
import os


# sudo chown thallock prefix -R

repos = [
	'billiards-attempts-api',
	'billiards-common',
	'billiards-graphics-api',
	'billiards-layouts-api',
	'billiards-scripts',
	'billiards-client',
	'billiards-config-api',
	'billiards-image-processing-api',
	'billiards-projection-api',
	'billiards-shots-api'
]


def get_repos_dir():
	if 'REPOS' in os.environ and os.path.exists(os.environ['REPOS']):
		repos_dir = os.environ['REPOS']
	else:
		repos_dir = os.path.join(__file__, '..', '..')
	return os.path.realpath(repos_dir)


def print_status(repo_dir):
	srcs = []
	repo = Repo(repo_dir)
	NUM_CHARS = 100
	print('#' + ''.join(['+'] * NUM_CHARS))
	print('#' + repo_dir)
	print('#' + ''.join(['+'] * NUM_CHARS))
	
	print('#\tmodified:')
	for item in repo.index.diff(None):
		if item.a_path.find('src/') == 0:
			srcs.append(item.a_path)
		print('#\t\t' + item.a_path)
	
	print('#\tuntracked:')
	for item in repo.untracked_files:
		if item.find('src/') == 0:
			srcs.append(item.a_path)
		print('#\t\t' + item)
	
	print('#\n#')
	return srcs


def main(message):
	repos_dir = get_repos_dir()
	cmds = []
	for directory in repos:
		repo_dir = os.path.join(repos_dir, directory)
		
		try:
			srcs = print_status(repo_dir)
			if len(srcs) == 0:
				continue
			
			cmds.append('pushd "' + repo_dir + "'")
			for src in srcs:
				cmds.append('git add "' + src + '"')
			cmds.append('git commit -m "' + message + '"')
			cmds.append('git push')
			cmds.append('popd')
		except:
			pass
	
	if message:
		for cmd in cmds:
			print(cmd)


if __name__ == '__main__':
	main(sys.argv[1] if len(sys.argv) >= 2 else None)
