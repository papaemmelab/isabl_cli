"""isabl_cli setup.py."""
from pathlib import Path
from setuptools import find_packages
from setuptools import setup


ROOT = Path(__file__).parent.resolve()
long_description = (ROOT / "README.md").read_text(encoding="utf-8")
version = (ROOT / "isabl_cli" / "VERSION").read_text(encoding="utf-8").strip()

setup(
    # single source package version
    version=version,
    # in combination with recursive-includes in MANIFEST.in, non-python files
    # within the isabl_cli will be copied into the
    # site-packages and wheels installation directories
    include_package_data=True,
    # return a list all Python packages found within the ROOT directory
    packages=find_packages(),
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Topic :: Utilities"
    ],
    entry_points={
        "console_scripts": ["isabl=isabl_cli.cli:main"]
    },
    setup_requires=["pytest-runner>=2.11.1"],
    install_requires=[
        "analytics-python>=1.4.0",
        "cached_property>=1.4.2",
        "Click>=7.0.0",
        "munch>=2.3.2",
        "python-slugify>=1.1.2",
        "pytz>=2018.4",
        "PyYAML>=5.3",
        "requests>=2.19.1",
        "six>=1.11.0",
        "packaging>=21.0"
    ],
    extras_require={
        "test": [
            "black>=18.9b0",
            "coverage>=4.4.2",
            "factory-boy>=2.11.1 ",
            "ipdb>=0.11",
            "pluggy>=0.7.1",
            "pydocstyle>=2.1.1",
            "pylint>=1.8.1",
            "pytest>=4.5",
            "pytest-cov>=2.5.1",
            "pytest-env>=0.6.2",
            "pytest-sugar>=0.9.1",
            "tox>=2.9.1"
        ]
    },
    author="Juan S. Medina",
    keywords=[],
    license="BSD",
    name="isabl-cli",
    test_suite="tests",
    long_description=long_description,
    long_description_content_type="text/markdown",
    description="isabl_cli command line client.",
    url="https://github.com/papaemmelab/isabl_cli"
)
