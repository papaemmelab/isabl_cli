# Isabl CLI

[![Run Isabl CLI tests][gh_actions_badge]][gh_actions_base]
[![pypi badge][pypi_badge]][pypi_base]
[![codecov badge][codecov_badge]][codecov_base]
[![code formatting][black_badge]][black_base]

> **_Terms of Use: Academic Research_**

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

## Usage

```bash
# Set up your environment variables
export ISABL_API_URL="http://<your-isabl-instance>/api/v1/"
export CLIENT_ID="<your-client-id>"

# Authenticate and cache credentials
isabl login

# See all available commands
isabl --help
```

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