# Isabl CLI

[![Run Isabl CLI tests][gh_actions_badge]][gh_actions_base]
[![pypi badge][pypi_badge]][pypi_base]
[![codecov badge][codecov_badge]][codecov_base]
[![code formatting][black_badge]][black_base]
[![MCP][mcp_badge]][mcp_base]

> **_Terms of Use: Academic Research_**

> [!TIP]
> Talk to your genomics data with AI. Try the [Isabl MCP Server](https://github.com/juanesarango/isabl-mcp) with Claude Code, Cursor, or your favorite agentic tool.

`isabl_cli` is a command line client and SDK for the Isabl Platform infrastructure.

## 🧬 Isabl Platform

- Isabl is a platform for the integration, management, and processing of individual-centric multimodal data. Welcome to the Isabl Documentation!

- Isabl is a plug-and-play data science framework designed to support the processing of multimodal patient-centric data. Have questions? Ask here.

- Isabl has been developed by the [Elli Papaemmanuil's Lab](https://papaemmelab.org) at [Memorial Sloan Kettering Cancer Center](https://www.mskcc.org/research/ski), and it's currently used by multiple scientific groups and contributors.

- Learn more: [Isabl documentation][documentation].

## 📓 Publication

>Medina-Martínez, J.S., Arango-Ossa, J.E., Levine, M.F. et al. Isabl Platform, a digital biobank for processing multimodal patient data. BMC Bioinformatics 21, 549 (2020).

<https://doi.org/10.1186/s12859-020-03879-7>

## Installation

```bash
pip install isabl_cli
```

## Development Setup

### Prerequisites

- Python 3.6 or higher
- Docker and Docker Compose (for running integration tests)
- Git (for cloning the API repository during tests)

### Setup Steps

1. Clone the repository:

   ```bash
   git clone https://github.com/papaemmelab/isabl_cli.git
   cd isabl_cli
   ```

2. Install the package in editable mode with test dependencies:

   ```bash
   pip install -e .[test]
   ```

   This installs the package in development mode along with all testing dependencies (pytest, coverage, black, etc.).

3. Verify the installation:

   ```bash
   isabl --help
   ```

## Testing

The test suite can be run in two ways:

### Running Tests via Script (Recommended)

The test script (`./tests/test-with-compose.sh`) automatically clones the `isabl_api` repository, starts the API and database through Docker containers, then executes the test suite using pytest with coverage reporting.

Run the full test suite:

```bash
./tests/test-with-compose.sh
```

To skip the build step (faster for subsequent runs):

```bash
./tests/test-with-compose.sh --skip-build
```

**Note**: The test script requires the API to be running. It will clone the `isabl_api` repository into a sibling `api` directory if it doesn't exist.

### Running Tests Directly

If you have a local Isabl API instance running (e.g., via `docker compose up -d` in the `isabl_api` repository), you can run pytest commands directly:

```bash
# Set required environment variables
export ISABL_CLIENT_ID=test-cli-client
export ISABL_API_URL=http://localhost:8000/api/v1/

# Run tests with coverage
pytest --cov=isabl_cli tests/
```

This allows for faster iteration during development without rebuilding containers.

### Running Specific Tests

You can run specific test files or test functions:

```bash
# Run a specific test file
pytest tests/test_api.py

# Run a specific test function
pytest tests/test_api.py::test_function_name

# Run with verbose output
pytest -v tests/
```

## Usage

```bash
# Set up your environment variables
export ISABL_API_URL="http://<your-isabl-instance>/api/v1/"
export CLIENT_ID="<your-client-id>"

# Authenticate and cache credentials
isabl login

# Single sign-on users: paste an API token instead of a password
# (copy it from the web UI user menu → "Copy API token")
isabl login --token <your-api-token>

# See all available commands
isabl --help
```

For scripted or CI environments, set `ISABL_API_TOKEN` to skip interactive
login entirely; the token is used directly and never written to disk.

## CI/CD Setup

The GitHub Actions workflow requires the following secrets to be configured in the repository settings:

### Required Secrets

- **`GH_PAT`**: GitHub Personal Access Token for authenticating with GitHub (used to clone the `isabl_api` repository during tests)
- **`DOCKER_USERNAME`**: Docker Hub username for the [papaemmelab organization](https://hub.docker.com/orgs/papaemmelab/repositories)
- **`DOCKER_PASSWORD`**: Docker Hub password or access token for the [papaemmelab organization](https://hub.docker.com/orgs/papaemmelab/repositories)
- **`CODECOV_TOKEN`**: Codecov token for the [papaemmelab/isabl_cli repository](https://app.codecov.io/gh/papaemmelab/isabl_cli)

Configure these secrets in your GitHub repository under **Settings → Secrets and variables → Actions**.

## Contributing

Contributions are welcome and greatly appreciated, check our [contributing guidelines]!

## License

isabl is free for academic and research use. For commercial use, please contact us.

## Credits

Isabl was created at [Elli Papaemmanuil's lab].

[documentation]: https://docs.isabl.io
[contributing guidelines]: https://docs.isabl.io/contributing-guide
[elli papaemmanuil's lab]: https://www.mskcc.org/research-areas/labs/elli-papaemmanuil
[black_badge]: https://img.shields.io/badge/code%20style-black-000000.svg
[black_base]: https://github.com/ambv/black
[codecov_badge]: https://codecov.io/github/papaemmelab/isabl_cli/graph/badge.svg?token=ODJ8DU73PH
[codecov_base]: https://codecov.io/github/papaemmelab/isabl_cli
[pypi_badge]: https://img.shields.io/pypi/v/isabl_cli.svg
[pypi_base]: https://pypi.python.org/pypi/isabl_cli
[gh_actions_badge]: https://github.com/papaemmelab/isabl_cli/actions/workflows/run-test.yaml/badge.svg
[gh_actions_base]: https://github.com/papaemmelab/isabl_cli/actions/workflows/run-test.yaml
[mcp_badge]: https://img.shields.io/badge/MCP-compatible-blue
[mcp_base]: https://github.com/juanesarango/isabl-mcp
