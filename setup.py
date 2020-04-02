import os
from setuptools import setup
from setuptools import find_packages


BINPREFIX = "nptt"


def gen_entry_points(pkgdir):
    """ Generate entry_points={'console_scripts': []} records
    relative to `pkgdir`.
    """
    results = []
    for root, dirs, files in os.walk(pkgdir):
        for f in files:
            if not f.endswith('.py') or f.endswith('__init__.py'):
                continue
            modname = os.path.splitext(f)[0]
            modpath = os.path.join(root, modname).replace(os.sep, ".")
            result = "{}_{}={}:{}".format(BINPREFIX, modname, modpath, "main")
            results += [result]

    return results


ENTRY_POINTS = gen_entry_points(os.path.normpath("nirspec_pipe_testing_tool/utils"))


setup(
    name="nirspec_pipe_testing_tool",
    use_scm_version=True,
    author="Maria Pena Guerrero",
    description="FILL THIS IN",
    url="https://github.com/spacetelescope/nirspec_pipe_testing_tool",
    license="BSD",
    classifiers=[
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Software Development :: Testing",
    ],
    python_requires=">=3.6",
    setup_requires=[
        "setuptools_scm",
    ],
    install_requires=[
        "astropy",
        "jwst @ git+https://github.com/spacetelescope/jwst#branch=master",
        "matplotlib",
        "numpy",
        "pytest",
    ],
    packages=find_packages(),
    entry_points = {
        "console_scripts": ENTRY_POINTS,
    },
)
