#!/bin/bash

# docker build -t billiards -f docker/Dockerfile docker


export REPOS=$(realpath ../)
export BUILD_DIR=$(realpath ./build/)
export PREFIX=$(realpath ./prefix/)

mkdir -p $BULD_DIR
mkdir -p $PREFIX


#  \
# 		billiards-graphics-api
		
for project in \
		billiards-layouts-api \
		billiards-projection-api \
		billiards-graphics-api \
		billiards-shots-api \
		billiards-config-api
do
	mkdir -p $BUILD_DIR/$project
	docker run --rm \
		-v $(realpath $REPOS/$project/):/source \
		-v $(realpath $BUILD_DIR/$project):/build \
		billiards \
		cmake /source -B /build

	docker run --rm \
		-v $(realpath $REPOS/$project/):/source \
		-v $(realpath $BUILD_DIR/$project):/build \
		-v $(realpath $REPOS/billiards-common/src/common):/usr/include/common \
		--workdir=/build \
		billiards \
		make

	docker run --rm \
		-v $(realpath $REPOS/$project/):/source \
		-v $(realpath $BUILD_DIR/$project):/build \
		-v $(realpath $REPOS/billiards-common/src/common):/usr/include/common \
		-v $(realpath prefix):/app \
		--workdir=/build \
		billiards \
		make DESTDIR=/app install
done


mkdir -p $BUILD_DIR/qt-projection
pushd $BUILD_DIR/qt-projection
DESTDIR=$PREFIX/native qmake $REPOS/billiards-projection-api
make
make install
popd

# docker run --rm \
# 	-v $(realpath $REPOS/billiards-projection-api/):/source \
# 	-v $(realpath $BUILD_DIR/qt-projection):/build \
# 	--env DESTDIR=$PREFIX \
# 	--workdir=/build \
# 	billiards \
# 	qmake-qt5 /source
# 
# docker run --rm \
# 	-v $(realpath $REPOS/billiards-projection-api/):/source \
# 	-v $(realpath $BUILD_DIR/qt-projection):/build \
# 	-v $(realpath $REPOS/billiards-common/src/common):/usr/include/common \
# 	--workdir=/build \
# 	billiards \
# 	make
# 
# docker run --rm \
# 	-v $(realpath $REPOS/billiards-projection-api/):/source \
# 	-v $(realpath $BUILD_DIR/qt-projection):/build \
# 	-v $(realpath $REPOS/billiards-common/src/common):/usr/include/common \
# 	-v $(realpath prefix):/app \
# 	--workdir=/build \
# 	billiards \
# 	make install
# 
# 
# docker run --rm \
# 	-v $(realpath prefix):/app \
# 	--workdir=/app \
# 	billiards \
# 	./qt-projection-api -platform xcb
