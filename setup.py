#!/usr/bin/env python3

from setuptools import setup, find_packages

# read the contents of your README file
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

packages = find_packages()
package_data = {"utils": ["utils/*"]}

data_files = [(".", ["LICENSE", "README.md"])]

setup(
    name="vstrains",
    version="1.1.0",
    zip_safe=True,
    author="Runpeng Luo and Yu Lin",
    author_email="runpengluo@gmail.com",
    description="VStrains: De Novo Reconstruction of Viral Strains via Iterative Path Extraction From Assembly Graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/metagentools/VStrains",
    license="MIT",
    packages=packages,
    package_data=package_data,
    data_files=data_files,
    include_package_data=True,
    scripts=["vstrains"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        # "graph-tool",
        "minimap2",
        "numpy",
        "gfapy",
        "matplotlib",
    ],
    python_requires=">=3",
)
