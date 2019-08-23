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
cd $API_DIR && docker-compose build && docker-compose up -d

# give some time to API to start and test Isabl CLI
echo "Giving 20 seconds to API to get started..."
sleep 20 && cd $CLI_DIR && pytest -v tests/
