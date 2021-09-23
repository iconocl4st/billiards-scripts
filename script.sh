#!/bin/bash

export REPOS=$(realpath ../..)
export PREFIX=$(realpath ../prefix)

for project in \
		billiards-common \
		billiards-layouts-api \
		billiards-projection-api \
		billiards-shots-api \
		billiards-config-api \
		billiards-graphics-api
do
	mkdir -p $project
	cmake $REPOS/$project -B $project
	pushd $project
	make
	make DESTDIR=$PREFIX install
	popd
done


mkdir -p qt-projection
pushd qt-projection
DESTDIR=$PREFIX qmake $REPOS/billiards-projection-api
make
make install
popd



# rm -rf qt-projection \
# 	billiards-layouts-api \
# 	billiards-projection-api \
# 	billiards-shots-api \
# 	billiards-config-api \
# 	billiards-graphics-api
