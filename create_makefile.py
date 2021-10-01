import sys
import os


# sudo chown thallock prefix -R

cmake_projects = [
	'billiards-layouts-api',
	'billiards-graphics-api',
	'billiards-shots-api',
	'billiards-config-api',
#	'billiards-projection-api',
]


clean_dirs = [
	'prefix/native',
	'prefix/app',
]


cleaner_dirs = [
	'build/billiards-graphics-api/',
	'build/billiards-layouts-api/',
	'build/billiards-projection-api/',
	'build/qt-projection/',
	'build/billiards-config-api/',
	'build/billiards-shots-api/'
]


def get_paths():
	if 'REPOS' in os.environ and os.path.exists(os.environ['REPOS']):
		repos_dir = os.environ['REPOS']
	else:
		repos_dir = os.path.join(__file__, '..', '..')
	repos_dir = os.path.realpath(repos_dir)
	scripts_dir = os.path.join(repos_dir, 'billiards-scripts')
	return {
		'repos': repos_dir,
		'scripts': scripts_dir,
		'prefix': os.path.join(scripts_dir, 'prefix'),
		'build': os.path.join(scripts_dir, 'build'),
		'makefile': os.path.join(scripts_dir, 'Makefile'),
		'common': os.path.join(repos_dir, 'billiards-common/src/common')
	}


def project_name(project):
	return project[len('billiards-'):-len('-api')]


def header_rule(keys):
	return '''
MAKEFLAGS += -j5
default: all

docker-image:
	docker build -t billiards -f docker/Dockerfile docker

.SILENT: run_cmd
run_cmd:
	echo "pushd {scripts} && docker-compose -f docker-compose.yml up"
	echo "pushd {scripts} && ./prefix/native/app/qt-projection-api"
	echo "pushd {repos}/billiards-client && npm start"

directories:
	mkdir -p "{prefix}"
	mkdir -p "{build}"

native-projection: directories
	mkdir -p "{build}/qt-projection"
	cd "{build}/qt-projection" && REPOS="{repos}" DESTDIR={prefix}/native qmake {repos}/billiards-projection-api
	cd "{build}/qt-projection" && make
	cd "{build}/qt-projection" && make install

update-makefile:
	python3 {scripts}/create_makefile.py
'''.format(**keys)


def project_rule(keys):
	return '''
.SILENT: {project-name}
{project-name}: directories docker-image
	mkdir -p {build}/{project}
	docker run --rm \\
		-v "{repos}/{project}/":/source \\
		-v "{build}/{project}/":/build \\
		billiards \\
		cmake /source -B /build

	docker run --rm \\
		-v "{repos}/{project}/":/source \\
		-v "{build}/{project}/":/build \\
		-v "{common}/":/usr/include/common \\
		--workdir=/build \\
		billiards \\
		make

	docker run --rm \\
		-v "{repos}/{project}/":/source \\
		-v "{build}/{project}/":/build \\
		-v "{common}/":/usr/include/common \\
		-v "{prefix}/":/app \\
		--workdir=/build \\
		billiards \\
		make DESTDIR=/app install

.SILENT: native-{project-name}
native-{project-name}: directories docker-image
	mkdir -p "{build}/native/{project}"
	cd "{repos}/{project}" && REPOS="{repos}" cmake -DCMAKE_BUILD_TYPE=Debug -B "{build}/native/{project}"
	cd "{build}/native/{project}" && make
	cd "{build}/native/{project}" && DESTDIR="{prefix}/native/" make install
'''.format(**keys)


def footer_rule(keys):
	return '''
clean:
	{clean-dirs}

cleaner: clean
	{cleaner-dirs}

all: native-projection {projects-string}
	'''.format(**keys)


def write_makefile(paths, outfile):
	outfile.write(header_rule(paths))
	for project in cmake_projects:
		outfile.write(project_rule({
			**paths,
			'project-name': project_name(project),
			'project': project
		}))
	outfile.write(footer_rule({
		**paths,
		'projects-string': ' '.join(project_name(p) for p in cmake_projects),
		'clean-dirs': '\n\t'.join('rm -rf "' + paths['scripts'] + '/' + p + '"' for p in clean_dirs),
		'cleaner-dirs': '\n\t'.join('rm -rf "' + paths['scripts'] + '/' + p + '"' for p in cleaner_dirs),
	}))


def main():
	paths = get_paths()
	with open(paths['makefile'], 'w') as outfile:
		write_makefile(paths, outfile)


if __name__ == '__main__':
	#MIN_PYTHON = (3, 1)
	#assert sys.version_info >= MIN_PYTHON, f"requires Python {'.'.join(str(n) for n in MIN_PYTHON)} or newer"
	main()
