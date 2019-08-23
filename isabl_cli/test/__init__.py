"""Isabl CLI testing utilities for apps development."""

import pytest


@pytest.fixture
def commit(request):
    """Return true if user wants to commit the test applications."""
    return request.config.getoption("--commit")


def pytest_addoption(parser):
    """Add option to commit or not applications."""
    parser.addoption("--commit", action="store_true", help="commit applications")
