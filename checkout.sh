
export REPOS=$(realpath ../../testing)

# git clone git@github.com:iconocl4st/billiards-scripts.git

pushd $REPOS
git clone git@github.com:iconocl4st/billiards-attempts-api.git
git clone git@github.com:iconocl4st/billiards-common.git
git clone git@github.com:iconocl4st/billiards-graphics-api.git
git clone git@github.com:iconocl4st/billiards-layouts-api.git
git clone git@github.com:iconocl4st/billiards-client.git
git clone git@github.com:iconocl4st/billiards-config-api.git
git clone git@github.com:iconocl4st/billiards-projection-api.git
git clone git@github.com:iconocl4st/billiards-shots-api.git
popd
