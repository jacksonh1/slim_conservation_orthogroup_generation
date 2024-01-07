import os

from setuptools import find_packages, setup

setup(
    name="local_orthogroup_tools",
    # mandatory
    version="0.1",
    # mandatory
    author="Jackson Halpin",
    packages=find_packages(
        where='src',
    ),
    package_dir={"": "src"}
)
