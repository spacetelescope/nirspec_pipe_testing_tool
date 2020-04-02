import os
from setuptools import setup
from setuptools import find_packages


def gen_entry_points(pkgdir, e_prefix):
    """ Generate entry_points={'console_scripts': []} records
    relative to `pkgdir`.
    """
    results = []
    root = pkgdir
    for f in os.listdir(pkgdir):
        # Skip sub-directories
        if os.path.isdir(f):
            continue
        # Skip non-script files and __init__
        if not f.endswith('.py') or f.endswith('__init__.py'):
            continue

        # Python module name is derived from the filename without its extension
        modname = os.path.splitext(f)[0]
        # Python module path is the absolute path using "." instead of "/"
        modpath = os.path.join(root, modname).replace(os.sep, ".")
        # Create record
        result = "{}_{}={}:{}".format(e_prefix, modname, modpath, "main")
        # Save record
        results += [result]

    return results


PACKAGE_NAME = "nirspec_pipe_testing_tool"
BINPREFIX = "nptt"
ENTRY_POINTS_PATH=os.path.normpath("{}/utils".format(PACKAGE_NAME))
ENTRY_POINTS = gen_entry_points(ENTRY_POINTS_PATH, BINPREFIX)


setup(
    name=PACKAGE_NAME,
    use_scm_version=True,
    author="Maria Pena Guerrero",
    description="FILL THIS IN",
    url="https://github.com/spacetelescope/{}".format(PACKAGE_NAME),
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
        "pysiaf",
        "pytest",
    ],
    packages=find_packages(),
    entry_points = {
        "console_scripts": ENTRY_POINTS,
    },
)
