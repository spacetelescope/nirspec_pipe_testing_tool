"""
py.test configuration for the *entire* test suite
"""

import pytest
import configparser


# HEADER
__author__ = "M. A. Pena-Guerrero"
__version__ = "1.0"

# HISTORY
# Nov 2017 - Version 1.0: initial version completed


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
    config = configparser.ConfigParser()
    config.read(request.config.getoption("--config_file"))
    # place the report in the working directory unless otherwise specified
    output_dir = config.get("calwebb_spec2_input_file", "output_directory")
    request.htmlpath = request.config.getoption('htmlpath', output_dir+"/report.html")
    config.read(request.config.getoption("--config_file"))
    return config


"""
@pytest.mark.hookwrapper
def pytest_runtest_makereport(item, call):
    pytest_html = item.config.pluginmanager.getplugin('html')
    outcome = yield
    report = outcome.get_result()
    extra = getattr(report, 'extra', [])
    if report.when == 'call':
        # get filename between square brackets
        m = re.match('^.*\[(.*)\].*$', item.name)
        fname = item.name.split('[')[0]+'_'+m.group(1).split('/')[-1][:-5]+'.pdf'
        # add plot if it exists
        if os.path.isfile(fname):
            extra.append(pytest_html.extras.image(fname))
        report.extra = extra
"""
