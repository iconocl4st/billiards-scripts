
from git import Repo
import sys
import os


# sudo chown thallock prefix -R

repos = [
	'billiards-attempts-api',
	'billiards-client',
	'billiards-common',
	'billiards-config-api',
	'billiards-graphics-api',
	'billiards-image-processing-api',
	'billiards-layouts-api',
	'billiards-projection-api',
	'billiards-scripts',
	'billiards-shots-api'
]


def get_repos_dir():
	if 'REPOS' in os.environ and os.path.exists(os.environ['REPOS']):
		repos_dir = os.environ['REPOS']
	else:
		repos_dir = os.path.join(__file__, '..', '..')
	return os.path.realpath(repos_dir)


def print_status(repo_dir, subdir):
	repo = Repo(repo_dir)
	
	modified = [item.a_path for item in repo.index.diff(None)]
	untracked = repo.untracked_files
	if len(modified) == 0 and len(untracked) == 0:
		return []
		
	NUM_CHARS = 100
	print('#' + ''.join(['+'] * NUM_CHARS))
	print('#' + repo_dir)
	print('#' + ''.join(['+'] * NUM_CHARS))
	
	
	print('#\tmodified:')
	for item in modified:
		print('#\t\t' + item)
	
	print('#\tuntracked:')
	for item in untracked:
		print('#\t\t' + item)
	
	print('#\n#')
	
	return [p for p in untracked + modified if p.find(subdir) == 0]


def main(message):
	repos_dir = get_repos_dir()
	cmds = []
	for directory in repos:
		repo_dir = os.path.join(repos_dir, directory)
		
		srcs = print_status(repo_dir, subdir='src/')
		if len(srcs) == 0:
			continue
		
		cmds.append('pushd "' + repo_dir + '"')
		for src in srcs:
			cmds.append('git add "' + src + '"')
		cmds.append('git commit -m "' + message + '"')
		cmds.append('git push')
		cmds.append('popd')
	
	if message:
		for cmd in cmds:
			print(cmd)


if __name__ == '__main__':
	# subdir
	main(sys.argv[1] if len(sys.argv) >= 2 else None)
