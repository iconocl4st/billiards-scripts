#!/bin/bash

for dire in $(ls):
do
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++"
	echo $dire
	echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++"
	pushd $dire
	git status
	popd
done
