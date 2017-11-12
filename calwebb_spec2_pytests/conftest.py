"""
py.test configuration for the *entire* test suite
"""

import pytest
import ConfigParser

def pytest_addoption(parser):
    """
    Specifies the files used for certain tests
    """
    parser.addoption("--config_file", action="store",
        help="specifies the file used for the test")
    parser.addoption("--gen_report", action="store_true",
        help="generate a report or not")

@pytest.fixture(scope="session", autouse=True)
def config(request):
    config = ConfigParser.ConfigParser()
    config.read(request.config.getoption("--config_file"))
    return config

