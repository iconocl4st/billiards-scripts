#!/bin/bash

export REPOS=$(realpath ../..)
ln -s $REPOS/billiards-common/ billiards-common
ln -s $REPOS/billiards-attempts-api/ billiards-attempts-api
ln -s $REPOS/billiards-graphics-api/ billiards-graphics-api
ln -s $REPOS/billiards-layouts-api/ billiards-layouts-api
ln -s $REPOS/billiards-projection-api/ billiards-projection-api
ln -s $REPOS/billiards-shots-api/ billiards-shots-api
ln -s $REPOS/billiards-config-api/ billiards-config-api


# rm common attempts graphics layouts projection shots
