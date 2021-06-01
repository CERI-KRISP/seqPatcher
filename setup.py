from sys import version
from setuptools import setup

setup(
    name="BackFill",
    version="0.1",
    scripts=["ab1_integration.py"],
    install_requires=[
        "pandas==1.2.4",
        "biopython==1.78",
        "click==7.1.2",
        "ruffus==2.8.4",
    ],
    # packages=[],
    # package_dir={},
    url="https://github.com/krisp-kwazulu-natal/sars-cov-2-sequencing-merge-sanger",
    license="GPLv3",
    author="Anmol Kiran; San James Emmanuel",
    author_email="anmol.kiran@gmail.com;sanemmanueljames@gmail.com",
    description="Inserts sars-cov2 S-gene sanger sequencing data in assemblies generated using HTS.",
)
