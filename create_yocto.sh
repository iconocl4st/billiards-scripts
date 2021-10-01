#!/bin/bash


mkdir -p yocto/output/dl_dir
mkdir -p yocto/output/sstate_dir
mkdir -p yocto/output/tmp_dir

mkdir -p yocto/sources
pushd yocto/sources
git clone git://git.yoctoproject.org/poky -b dunfell
git clone git://git.yoctoproject.org/meta-raspberrypi -b dunfell
git clone https://github.com/meta-qt5/meta-qt5 -b dunfell
git clone git://github.com/jumpnow/meta-rpi -b dunfell 
git clone https://git.openembedded.org/meta-openembedded -b dunfell
git clone https://github.com/jumpnow/meta-jumpnow -b dunfell
git clone git://git.yoctoproject.org/meta-security.git -b dunfell
git clone git@github.com:iconocl4st/meta-billiards.git
popd


mkdir -p yocto/build
cp -r yocto/conf yocto/build/
