# setup.py
import os
from setuptools import setup, find_packages

# Load version from version.py
version = {}
with open(os.path.join('vntyper', 'version.py')) as f:
    exec(f.read(), version)

setup(
    name="vntyper",
    version=version['__version__'],
    packages=find_packages(),
    include_package_data=True,  # Include package data as specified in MANIFEST.in or package_data
    install_requires=[
        "pandas>=2.2.0",
        "numpy>=2.0.2",
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
    author="Bernt Popp",
    author_email="bernt.popp.md@gmail.com",
    description="VNtyper: A tool for genotyping MUC1-VNTR",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/berntpopp/vntyper",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.9',
    extras_require={
        "dev": [
            "pytest",
            "black",
            "flake8",
        ],
    },
    package_data={
        'vntyper.scripts': ['kestrel_filter_config.json'],
        'vntyper': [
            'config.json',  # Include config.json in the vntyper package
            'templates/report_template.html',  # Include report_template.html
            'templates/cohort_summary_template.html',  # Include cohort_summary_template.html
        ],
    },
)
