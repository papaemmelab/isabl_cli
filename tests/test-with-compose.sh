# parse arguments
SKIP_BUILD=false
for arg in "$@"; do
  case $arg in
    --skip-build)
    SKIP_BUILD=true
    shift
    ;;
  esac
done

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
<<<<<<< HEAD
GH_BRANCH=${GITHUB_REF#refs/heads/}
=======
GH_BRANCH=$(echo ${GITHUB_REF#refs/heads/})
>>>>>>> 83fff43 (🔐 try appending GH_PAT to repo url)
ISABL_BRANCH=${GH_BRANCH:-master}
echo "ISABL_BRANCH set to $ISABL_BRANCH given Github branch: $GH_BRANCH"

<<<<<<< HEAD
# reset environment variables
unset ISABL_CLIENT_ID
unset ISABL_API_URL
=======
# clone api from github
rm -rf $API_DIR && git clone https://${GH_PAT}@github.com/papaemmelab/isabl_api.git $API_DIR
cd $API_DIR && (git checkout $ISABL_BRANCH || true) && docker-compose build && docker-compose up -d

# give some time to API to start and test Isabl CLI
echo "Giving 30 seconds to API to get started..."
>>>>>>> 8cb582f (🔐 try appending GH_PAT to repo url)
export ISABL_CLIENT_ID=test-cli-client

if [ "$SKIP_BUILD" = false ]; then
   # clone api from github
   rm -rf $API_DIR && git clone https://github.com/papaemmelab/isabl_api.git $API_DIR
   cd $API_DIR && (git checkout $ISABL_BRANCH || true) && docker-compose build && docker-compose up -d

   # give some time to API to start and test Isabl CLI
   echo "Giving 30 seconds to API to get started..."
   sleep 30
fi

cd $CLI_DIR && pytest -vs tests/ --cov=isabl_cli || (
   # print django logs if tests fail
   cd $API_DIR && docker-compose logs django && exit 1
)
