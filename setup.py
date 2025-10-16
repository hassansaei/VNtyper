# setup.py

import os
from setuptools import setup, find_packages

# Load version from version.py
version = {}
with open(os.path.join("vntyper", "version.py")) as f:
    exec(f.read(), version)

# ===========================
# Installation Reference
# ===========================
# After installing vntyper, you need to download and set up the required reference files.
# Run the following command to perform the installation of references:
#
#     vntyper --config-path /path/to/config.json install-references --output-dir /path/to/install
#
# For detailed instructions, refer to the README.md or visit:
# https://github.com/hassansaei/vntyper#installation

setup(
    name="vntyper",
    version=version["__version__"],
    packages=find_packages(),
    include_package_data=True,  # Include package data as specified in MANIFEST.in or package_data
    install_requires=[
        "pandas>=2.2.0,<2.3.0",
        "numpy>=1.26.0,<2.0.0",
        "regex>=2024.7.24",
        "biopython>=1.84",
        "setuptools>=72.2.0",
        "pysam>=0.22.1",
    ],
    entry_points={
        "console_scripts": [
            "vntyper=vntyper.cli:main",
        ],
    },
    author="Hassan Saei, Bernt Popp",
    author_email="hassan.saei@inserm.fr, bernt.popp.md@gmail.com",
    description="VNtyper: A tool for genotyping MUC1-VNTR",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/hassansaei/vntyper",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    extras_require={
        "dev": [
            "pytest",
            "black",
            "flake8",
        ],
    },
    package_data={
        "vntyper.scripts": [
            "kestrel_config.json",
            "install_references_config.json",
            "report_config.json",
        ],  # Include kestrel_config.json, install_references_config.json, and report_config.json
        "vntyper": [
            "config.json",  # Include config.json in the vntyper package
            "templates/report_template.html",  # Include report_template.html
            "templates/cohort_summary_template.html",  # Include cohort_summary_template.html
            "dependencies/kestrel/*.jar",  # Include all JAR files in dependencies/kestrel/
            "modules/shark/shark_config.json",  # Include shark_config.json
            "modules/advntr/advntr_config.json",  # Include advntr_config.json
        ],
    },
)
