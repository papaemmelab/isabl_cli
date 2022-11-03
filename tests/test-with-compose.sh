# get path to demo directory
if [ "$SHELL" = "zsh" ]; then
   TEST_DIR="$( cd "$(dirname ${(%):-%N})" ; pwd -P )"
else
   TEST_DIR="$( cd "$(dirname "$BASH_SOURCE")" ; pwd -P )";
fi

# set directories
CLI_DIR=${TEST_DIR}/..
API_DIR=${CLI_DIR}/api
echo "API directory set to: $API_DIR"
echo "CLI directory set to: $CLI_DIR"

# get current isabl branch in case this test depends on a particular branch
ISABL_BRANCH=${TRAVIS_PULL_REQUEST_BRANCH:-${TRAVIS_BRANCH:-master}}
echo "ISABL_BRANCH set to $ISABL_BRANCH given travis branch: $TRAVIS_BRANCH $TRAVIS_PULL_REQUEST_BRANCH"

# clone api from github
rm -rf $API_DIR && git clone git@github.com:papaemmelab/isabl_api.git $API_DIR
cd $API_DIR && (git checkout $ISABL_BRANCH || true) && docker-compose build && docker-compose up -d

# give some time to API to start and test Isabl CLI
echo "Giving 30 seconds to API to get started..."
export ISABL_CLIENT_ID=test-cli-client
sleep 30 && cd $CLI_DIR && pytest -vs tests/ --cov=isabl_cli || (
   # print django logs if tests fail
   cd $API_DIR && docker-compose logs django && exit 1
)
