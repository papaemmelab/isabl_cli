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

# clone api from github
[ ! -d $API_DIR ] && git clone git@github.com:isabl-io/api.git $API_DIR

# build container
[[ -z "${TRAVIS_BRANCH}" ]] && API_BRANCH='master' || API_BRANCH="${TRAVIS_BRANCH}"
cd $API_DIR && (git checkout $API_BRANCH || true) && docker-compose build && docker-compose up -d

# give some time to API to start and test Isabl CLI
echo "Giving 30 seconds to API to get started..."
sleep 30 && cd $CLI_DIR && pytest -v tests/ --cov=isabl_cli
