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

# create superuser
docker-compose run --rm django python manage.py shell -c "from django.contrib.auth.models import User; User.objects.create_superuser('admin', 'admin@admin.admin', 'admin')"

# test Isabl CLI
cd $CLI_DIR && pytest -v tests/
