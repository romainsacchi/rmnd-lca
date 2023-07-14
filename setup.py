import os
from pathlib import Path

from setuptools import setup

packages = []
root_dir = os.path.dirname(__file__)
if root_dir:
    os.chdir(root_dir)

# read the contents of your README file
this_directory = Path(__file__).parent
README = (this_directory / "README.md").read_text()

# Probably should be changed, __init__.py is no longer required for Python 3
for dirpath, dirnames, filenames in os.walk("premise"):
    # Ignore dirnames that start with '.'
    if "__init__.py" in filenames:
        pkg = dirpath.replace(os.path.sep, ".")
        if os.path.altsep:
            pkg = pkg.replace(os.path.altsep, ".")
        packages.append(pkg)


def package_files(directory):
    paths = []
    for path, directories, filenames in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join("..", path, filename))
    return paths


setup(
    name="premise",
    version="1.5.8",
    python_requires=">=3.9,<3.11",
    packages=packages,
    author="Romain Sacchi <romain.sacchi@psi.ch>, Alois Dirnaichner <dirnaichner@pik-potsdam.de>, Chris Mutel "
    "<chris.mutel@psi.ch>",
    license=open("LICENSE").read(),
    # Only if you have non-python data (CSV, etc.). Might need to change the directory name as well.
    include_package_data=True,
    install_requires=[
        "numpy",
        "wurst==0.3.3",
        "bw2io==0.8.7",
        "pandas==1.5.3",
        "bw2data",
        "xarray",
        "prettytable",
        "pycountry",
        "cryptography",
        "premise_gwp",
        "pyYaml",
        "sparse",
        "schema",
        "datapackage",
        "requests",
        "bottleneck",
        "constructive_geometries>=0.8.2",
    ],
    url="https://github.com/polca/premise",
    description="Coupling IAM output to ecoinvent LCA database ecoinvent for prospective LCA",
    long_description_content_type="text/markdown",
    long_description=README,
    classifiers=[
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
