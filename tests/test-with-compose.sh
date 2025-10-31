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

# require GitHub token (prefer GH_PAT for private repos, fallback to GITHUB_TOKEN)
if [ -n "$GH_PAT" ]; then
  TOKEN="$GH_PAT"
  TOKEN_NAME="GH_PAT"
  echo "Using GH_PAT for authentication (recommended for private repos)"
elif [ -n "$GITHUB_TOKEN" ]; then
  TOKEN="$GITHUB_TOKEN"
  TOKEN_NAME="GITHUB_TOKEN"
  echo "Using GITHUB_TOKEN for authentication"
  echo "NOTE: If cloning a private repo, GITHUB_TOKEN may not work. Use GH_PAT instead."
else
  echo "ERROR: Either GH_PAT or GITHUB_TOKEN is required for cloning dependencies." >&2
  echo "  - GH_PAT: Use this for private repositories (Personal Access Token)" >&2
  echo "  - GITHUB_TOKEN: Auto-provided by GitHub Actions (limited to same org)" >&2
  exit 1
fi

# debug: verify token is set (without exposing it)
echo "${TOKEN_NAME} is set (length: ${#TOKEN} characters)"

# configure netrc for git authentication
umask 077
printf "machine github.com\n  login x-access-token\n  password %s\n" "$TOKEN" > ~/.netrc

# debug: verify netrc was created
if [ -f ~/.netrc ]; then
  echo ".netrc file created successfully (permissions: $(stat -c '%a' ~/.netrc 2>/dev/null || stat -f '%OLp' ~/.netrc))"
  # show first line only (without password) for verification
  head -1 ~/.netrc | sed 's/password.*/password ***/'
else
  echo "ERROR: Failed to create .netrc file" >&2
  exit 1
fi

# configure git to use netrc
export GIT_TERMINAL_PROMPT=0
export GIT_ASKPASS=echo

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
GH_BRANCH=${GITHUB_REF#refs/heads/}
ISABL_BRANCH=${GH_BRANCH:-master}
echo "ISABL_BRANCH set to $ISABL_BRANCH given Github branch: $GH_BRANCH"

# reset environment variables
unset ISABL_CLIENT_ID
unset ISABL_API_URL
export ISABL_CLIENT_ID=test-cli-client

if [ "$SKIP_BUILD" = false ]; then
   # clone api from github
   echo "Attempting to clone papaemmelab/isabl_api..."
   echo "Target directory: $API_DIR"
   echo "Branch: $ISABL_BRANCH"
   rm -rf $API_DIR
   
   # test repo access first with a dry-run
   echo "Testing repository access..."
   if git ls-remote https://github.com/papaemmelab/isabl_api.git > /dev/null 2>&1; then
     echo "Repository access successful"
   else
     echo "ERROR: Cannot access repository. This could mean:" >&2
     echo "  1. The repository is private and GITHUB_TOKEN doesn't have permission" >&2
     echo "  2. The repository doesn't exist" >&2
     echo "  3. Authentication failed" >&2
     echo "" >&2
     echo "Debug info:" >&2
     echo "  Token used: ${TOKEN_NAME}" >&2
     echo "  Token length: ${#TOKEN}" >&2
     echo "  .netrc exists: $([ -f ~/.netrc ] && echo yes || echo no)" >&2
     echo "" >&2
     echo "If the repo is PRIVATE, you likely need to:" >&2
     echo "  1. Use GH_PAT (Personal Access Token) instead of GITHUB_TOKEN" >&2
     echo "  2. Set GH_PAT secret in GitHub repository settings" >&2
     git ls-remote https://github.com/papaemmelab/isabl_api.git 2>&1
     exit 1
   fi
   
   git clone https://github.com/papaemmelab/isabl_api.git $API_DIR
   if [ $? -ne 0 ]; then
     echo "ERROR: git clone failed" >&2
     exit 1
   fi
   cd $API_DIR && (git checkout $ISABL_BRANCH || true) && docker-compose build && docker-compose up -d

   # give some time to API to start and test Isabl CLI
   echo "Giving 30 seconds to API to get started..."
   sleep 30
fi

cd $CLI_DIR && pytest -vs tests/ --cov=isabl_cli || (
   # print django logs if tests fail
   cd $API_DIR && docker-compose logs django && exit 1
)
